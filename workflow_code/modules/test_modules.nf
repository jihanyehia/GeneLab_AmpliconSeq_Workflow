/*
 * Test-mode modules for AmpIllumina pipeline
 * 
 * Sample selection via composite quality score ensures representative,
 * non-outlier samples flow through all downstream test steps.
 *
 * Three levels (controlled by --test_level):
 *   primers          — select samples → detect primers → cutadapt grid (if primers found); stop
 *   filter           — select samples → filterAndTrim grid; stop (skips primer check)
 *   full (default)   — select samples → detect primers →
 *                        if found: cutadapt grid → filterAndTrim grid
 *                        if not:   filterAndTrim grid directly
 */

// ───────────────────────────────────────────────────────────────────────────────────────
// 1. SELECT SAMPLES - only if --input_file is not provided
//    Reads multiqc_fastqc.txt from the MultiQC data directory to compute a
//    composite score per sample from multiple metrics (mean quality, read count,
//    GC distance from median, duplication rate, adapter content). 
//    Ranks samples on the composite and selects N at evenly-spaced interior quantiles,
//    avoiding outliers at either extreme.
//    Emits:
//      selected_ids  — one selected sample ID per line (used to filter channels)
//      report        — full ranking TSV included in the HTML report
// ───────────────────────────────────────────────────────────────────────────────────────
process TEST_SELECT_SAMPLES {

    input:
    path multiqc_data_dir   // multiqc_data/ directory from RAW_MULTIQC
    val  n_samples
    val  assay_suffix

    output:
    path "selected_ids.txt",        emit: selected_ids
    path "sample_quality_rank.tsv", emit: report

    script:
    """
    test_select_samples.py \\
        --multiqc-data "${multiqc_data_dir}" \\
        --n            "${n_samples}" \\
        --assay-suffix "${assay_suffix}" \\
        --output       "sample_quality_rank.tsv" \\
        > selected_ids.txt
    """
}


// ───────────────────────────────────────────────────────────────────────────────────────
// 2a. PRIMER DETECTION
//    Searches the first N reads of each selected sample for the F/R primers
//    and their reverse complements 
//    Emits a per-sample TSV of hit rates
// ───────────────────────────────────────────────────────────────────────────────────────
process TEST_PRIMER_DETECTION {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads), val(isPaired)
    val   F_primer
    val   R_primer
    val   n_reads

    output:
    path "${sample_id}_primer_hits.tsv", emit: hits

    script:
    def r1 = reads instanceof List ? reads[0] : reads
    def r2 = reads instanceof List && reads.size() > 1 ? reads[1] : null
    def r2_arg = r2 ? "--r2 \"${r2}\" \\" : ""
    """
    test_primer_detection.py \\
        --sample-id "${sample_id}" \\
        --r1 "${r1}" \\
        ${r2_arg}
        --f-primer "${F_primer}" \\
        --r-primer "${R_primer}" \\
        --n-reads ${n_reads}
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// 2b. PARSE PRIMER RESULTS
//    Aggregates per-sample primer hit TSVs and emits a single-line decision
//    file ("primers_found" or "primers_not_found") used to branch the workflow.
// ─────────────────────────────────────────────────────────────────────────────
process TEST_PARSE_PRIMER_RESULTS {

    input:
    path primer_hits   // collected TSVs from TEST_PRIMER_DETECTION
    val  threshold     // min mean F_primer_pct to call primers present

    output:
    path "primer_detection_summary.tsv", emit: summary
    path "primer_decision.txt",          emit: decision

    script:
    """
    test_parse_primer_results.py \\
        --threshold       "${threshold}" \\
        --output-summary  "primer_detection_summary.tsv" \\
        --output-decision "primer_decision.txt"
    """
}


// ───────────────────────────────────────────────────────────────────────────────────────────────
// 3a. CUTADAPT PARAMETER GRID
//    Runs cutadapt with all combinations of anchored/unanchored × linked/unlinked
//    × discard_untrimmed on each selected sample, then parses the JSON stats into a one-row TSV.
// ───────────────────────────────────────────────────────────────────────────────────────────────
process TEST_CUTADAPT_GRID {

    tag "$sample_id — ${combo[0]}"

    input:
    tuple val(sample_id), path(reads), val(isPaired)
    val   F_primer
    val   R_primer
    val   min_len
    each  combo                  // [label, anchored_flag, linked_flag, discard_untrimmed_flag]

    output:
    tuple val(sample_id), path("*_trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_${combo[0]}_cutadapt_grid.tsv", emit: grid_row

    script:
    def combo_label   = combo[0]
    def anchored      = combo[1]  // true/false
    def linked        = combo[2]
    def discard       = combo[3]

    def r1 = reads instanceof List ? reads[0] : reads
    def r2 = isPaired == "true" ? reads[1] : ""

    def discardFlag = discard ? "--discard-untrimmed" : ""
    def json_out     = "${sample_id}_${combo_label}_cutadapt.json"
    """
    # Build adapter strings    
    F_RC=`echo ${F_primer} |tr ATGCRYSWKMBVDHN  TACGYRSWMKVBHDN |rev`
    R_RC=`echo ${R_primer} |tr ATGCRYSWKMBVDHN  TACGYRSWMKVBHDN |rev`
    
    # Anchored: ^PRIMER  Unanchored: PRIMER
    ANCHOR=""
    if [ "${anchored}" = "true" ]; then
        ANCHOR="^"
    fi

    # Linked (PE): -a PRIMER...RC(reverse complement of other primer)  Unlinked: -g PRIMER
    if [ "${linked}" = "true" ]; then
        # Linked adapters: forward adapter...RC(reverse adapter) on R1, vice versa on R2
        ADAPTER_R1="-a \${ANCHOR}${F_primer}...\${R_RC}"
        ADAPTER_R2="-A \${ANCHOR}${R_primer}...\${F_RC}"
    else
        ADAPTER_R1="-g \${ANCHOR}${F_primer}"
        ADAPTER_R2="-G \${ANCHOR}${R_primer}"
    fi

    if [ ${isPaired} = "true" ]; then
        cutadapt \\
            \$ADAPTER_R1 \$ADAPTER_R2 \\
            -m ${min_len} \\
            ${discardFlag} \\
            --json="${json_out}" \\
            -o "${sample_id}_${combo_label}_R1_trimmed.fastq.gz" \\
            -p "${sample_id}_${combo_label}_R2_trimmed.fastq.gz" \\
            "${r1}" "${r2}" \
            > /dev/null 2>&1 || true
    else
        cutadapt \\
            \$ADAPTER_R1 \\
            -m ${min_len} \\
            ${discardFlag} \\
            --json="${json_out}" \\
            -o "${sample_id}_${combo_label}_trimmed.fastq.gz" \\
            "${r1}" \\
            > /dev/null 2>&1 || true
    fi

    test_cutadapt_grid.py \\
        --sample-id "${sample_id}" \\
        --combo-label "${combo_label}" \\
        --anchored "${anchored}" \\
        --linked "${linked}" \\
        --discard "${discard}" \\
        --json-file "${json_out}"
    """
}


// ───────────────────────────────────────────────────────────────────────────────────────────────
// 3b. PICK BEST CUTADAPT COMBO
//     Parses all cutadapt grid TSVs to pick the best combo for each sample, then
//     emits a summary TSV and two text files with the best combo for discard vs keep.
//     If no cutadapt grid results are found (e.g. because no primers were detected),
//     emits placeholder files to avoid downstream errors.
// ───────────────────────────────────────────────────────────────────────────────────────────────

process TEST_PICK_BEST_CUTADAPT {

    input:
    path grid_rows   // all *_cutadapt_grid.tsv files collected from TEST_CUTADAPT_GRID

    output:
    path "best_cutadapt_combos.tsv", emit: summary
    path "best_discard_combo.txt",   emit: best_discard
    path "best_keep_combo.txt",      emit: best_keep
    path "raw_read_counts.tsv",      emit: raw_read_counts

    script:
    """
    test_pick_best_cutadapt.py \\
        --output-summary      "best_cutadapt_combos.tsv" \\
        --output-best-discard "best_discard_combo.txt" \\
        --output-best-keep    "best_keep_combo.txt"

    # Extract unique per-sample raw read counts from the cutadapt grid TSVs.
    # Column layout: sample, combo, anchored, linked, discard_untrimmed, reads_in, ...
    { printf 'sample\\treads_in\\n'; \\
      awk -F'\\t' 'FNR>1 && !seen[\$1]++{print \$1"\\t"\$6}' *_cutadapt_grid.tsv | sort; \\
    } > raw_read_counts.tsv
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// 4a. PICK TRUNCLEN
//     Uses the MultiQC JSON to pick a truncLen for each read direction based on
//     from raw MultiQC per-base sequence quality data. 
//     Emits a text file with the recommended truncLen for each read direction.
// ─────────────────────────────────────────────────────────────────────────────
process TEST_PICK_TRUNCLEN {

    input:
    path(multiqc_dir)
    val(isPaired)

    output:
    path "auto_trunclen.txt"

    script:
    """
    test_pick_truncLen.py \\
        --multiqc-json "${multiqc_dir}/multiqc_data.json" \\
        --is-paired    "${isPaired}" \\
        --output       "auto_trunclen.txt"
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// 4b. FILTER PARAMETER GRID
//    Runs DADA2::filterAndTrim for one truncLen / maxEE combination on the
//    selected sample reads. (One R process per parameter combination)
// ─────────────────────────────────────────────────────────────────────────────
process TEST_FILTER_GRID {

    tag "trunc_L${combo[0]}_R${combo[1]}__maxEE_L${combo[2]}_R${combo[3]}__${label}"

    input:
    tuple path(sample_IDs_file), val(isPaired)
    // label: trim status (e.g. best_discard, best_keep, or "raw"), reads: flat list of fastq.gz files staged into work dir by Nextflow, combo: [left_trunc, right_trunc, left_maxEE, right_maxEE, truncLen_source]
    tuple val(label), path(reads), val(combo)
    path raw_read_counts_file   // raw_read_counts.tsv from TEST_PICK_BEST_CUTADAPT, or NO_FILE when cutadapt was not run
    
    output:
    path "filter_grid_${combo[0]}_${combo[1]}_${combo[2]}_${combo[3]}_${label}_read_counts_tracking.tsv", emit: read_counts_tracking

    script:
    def (lt, rt, lm, rm, source) = combo
    """
    test_filter_grid.R \\
        "${params.assay_suffix}" \\
        "${sample_IDs_file}" \\
        "${isPaired}" \\
        "${label}" \\
        "${lt}" \\
        "${rt}" \\
        "${lm}" \\
        "${rm}" \\
        "${source}" \\
        "${params.concatenate_reads_only}" \\
        "filter_grid_${combo[0]}_${combo[1]}_${combo[2]}_${combo[3]}_${label}_read_counts_tracking.tsv"
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// 5. TEST REPORT
//    Collects all TSVs and renders a single self-contained HTML report
//    Publishes to outdir/Test_Mode_Report/
// ─────────────────────────────────────────────────────────────────────────────
process TEST_REPORT {

    input:
    path quality_rank    // sample_quality_rank.tsv from TEST_SELECT_SAMPLES
    path primer_summary  // primer_detection_summary.tsv from TEST_PARSE_PRIMER_RESULTS (or NO_FILE)
    path primer_decision // primer_decision.txt from TEST_PARSE_PRIMER_RESULTS
    path cutadapt_grid   // TSV(s) from TEST_CUTADAPT_GRID  (or NO_FILE)
    path best_cutadapt_combos // TSV from TEST_PICK_BEST_CUTADAPT (or NO_FILE)
    path filter_grid    // TSV(s) from TEST_FILTER_GRID    (or NO_FILE)
    val  test_level
    val  F_primer
    val  R_primer
    val  target_region
    val  accession_or_file
    path multiqc_dir
    val  is_paired
    val  auto_lt
    val  auto_rt

    output:
    path "${accession_or_file}_${target_region}_test_mode_report.html", emit: report

    script:
    def trunc_args = auto_lt == 'none' ? '' : "--auto-lt ${auto_lt} --auto-rt ${auto_rt}"
    """
    test_report.py \\
        --test-level "${test_level}" \\
        --f-primer "${F_primer}" \\
        --r-primer "${R_primer}" \\
        --target-region "${target_region}" \\
        --input "${accession_or_file}" \\
        --multiqc-json "${multiqc_dir}/multiqc_data.json" \\
        --is-paired "${is_paired}" \\
        ${trunc_args}
    """
}