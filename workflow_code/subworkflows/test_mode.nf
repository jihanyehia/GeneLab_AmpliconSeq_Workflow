/*
 * TEST_MODE subworkflow
 *
 * Called from main.nf when params.test_mode == true.
 * Staging and raw QC (including MultiQC) run on the full dataset.
 * If params.input_file is provided, this subworkflow uses that sample set directly.
 * Otherwise it selects representative samples from the full dataset using MultiQC.
 *
 * Flow:
 *   TEST_PRIMER_DETECTION -> TEST_PARSE_PRIMER_RESULTS
 *     │
 *     ├─ primers_found:     TEST_CUTADAPT_GRID -> TEST_PICK_BEST_CUTADAPT
 *     │                         │
 *     │                         └─ TEST_PICK_TRUNCLEN -> TEST_FILTER_GRID -> TEST_REPORT
 *     │
 *     └─ primers_not_found: TEST_PICK_TRUNCLEN -> TEST_FILTER_GRID directly -> TEST_REPORT
 *
 * Parameters consumed:
 *   params.test_level            — 'primers' | 'filter' | 'full'  (default: 'full')
 *   params.test_n_samples        — number of representative samples to select when input_file is not provided
 *   params.test_primer_nreads    — reads to inspect for primer detection (default: 10000)
 *   params.test_primer_threshold — min mean F_primer_pct to call primers present (default: 50.0)
 *   params.test_trunc_left       — optional comma-sep list of extra left truncLen values (default: null, auto-pick only)
 *   params.test_trunc_right      — optional comma-sep list of extra right truncLen values (default: null, auto-pick only)
 *   params.test_maxEE_left       — comma-sep list of left maxEE values  e.g. "1,2"
 *   params.test_maxEE_right      — comma-sep list of right maxEE values e.g. "1,2"
 *
 * Outputs:
 *   report  — path to test_mode_report.html (published to outdir/Test_Mode_Report/)
 */

include { TEST_SELECT_SAMPLES       } from '../modules/test_modules.nf'
include { TEST_PRIMER_DETECTION     } from '../modules/test_modules.nf'
include { TEST_PARSE_PRIMER_RESULTS } from '../modules/test_modules.nf'
include { TEST_CUTADAPT_GRID        } from '../modules/test_modules.nf'
include { TEST_PICK_BEST_CUTADAPT   } from '../modules/test_modules.nf'
include { TEST_PICK_TRUNCLEN        } from '../modules/test_modules.nf'
include { TEST_FILTER_GRID          } from '../modules/test_modules.nf'
include { TEST_REPORT               } from '../modules/test_modules.nf'

workflow TEST_MODE {

    take:
    staged_reads_ch   // tuple(sample_id, [reads], isPaired) from STAGE_ANALYSIS (full dataset or samples specified in --input_file)
    runsheet_ch       // path to runsheet.csv from STAGE_ANALYSIS (full dataset or samples specified in --input_file)
    multiqc_data_ch   // path to multiqc_data/ directory from RAW_MULTIQC
    primers_ch        // tuple(F_primer, R_primer) from STAGE_ANALYSIS

    main:

    def test_level = params.test_level ?: 'full'
    def n_reads        = params.test_primer_nreads    ?: 10000
    def threshold      = params.test_primer_threshold ?: 5.0
    
    // ─────────────────────────────────────────────────────────────────────────
    // Step 1 — Sample selection (all levels)
    // If --input_file is provided, use that sample set directly
    // ─────────────────────────────────────────────────────────────────────────    
    if (!params.input_file) {
        
        TEST_SELECT_SAMPLES(multiqc_data_ch, params.test_n_samples ?: 3, params.assay_suffix ?: '')

        // Emits flat string IDs
        selected_ids_ch = TEST_SELECT_SAMPLES.out.selected_ids
            .splitText()
            .map { it.trim() }
            .filter { it }

        // Use standard cross-referencing operator
        selected_reads_ch = staged_reads_ch
            .join(selected_ids_ch) // Matches on sample_id by default
            .map { sample_id, reads, isPaired -> tuple(sample_id, reads, isPaired) }

        sample_quality_rank_ch = TEST_SELECT_SAMPLES.out.report
    }
    else {
        // If --input_file is provided, bypass sample selection and use the provided sample set directly
        selected_reads_ch = staged_reads_ch
        sample_quality_rank_ch = channel.of(file('NO_FILE_QUALITY_RANK'))
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Step 2 — primer detection (primers and full levels only)
    // ─────────────────────────────────────────────────────────────────────────
    F_primer_ch = primers_ch.map { it[0] }
    R_primer_ch = primers_ch.map { it[1] }
    if (test_level == 'primers' || test_level == 'full'){

        TEST_PRIMER_DETECTION(
            selected_reads_ch,
            F_primer_ch,
            R_primer_ch,
            n_reads
        )

        TEST_PARSE_PRIMER_RESULTS(
            TEST_PRIMER_DETECTION.out.hits.collect(),
            threshold
        )
        primer_summary_ch = TEST_PARSE_PRIMER_RESULTS.out.summary
        primer_decision_ch = TEST_PARSE_PRIMER_RESULTS.out.decision

        
        // Only run cutadapt grid when primers were detected in the sample set
        // Otherwise skip to report with a note that primers were not found
        primers_found_reads_ch = selected_reads_ch
                .combine(TEST_PARSE_PRIMER_RESULTS.out.decision)
                .filter { sample_id, reads, isPaired, decision_file ->
                    decision_file.text.trim() == 'primers_found'
                }
                .map { sample_id, reads, isPaired, decision_file ->
                    tuple(sample_id, reads, isPaired)
                }
        
        // ── cutadapt combos definition ────────────────────────────────────────
        // [label, anchored, linked, discard_untrimmed]
        def cutadapt_combos = [
            ['anchored_linked_discard',    true,  true,  true ],
            ['anchored_linked_keep',       true,  true,  false],
            ['anchored_unlinked_discard',  true,  false, true ],
            ['anchored_unlinked_keep',     true,  false, false],
            ['unanchored_linked_discard',  false, true,  true ],
            ['unanchored_linked_keep',     false, true,  false],
            ['unanchored_unlinked_discard',false, false, true ],
            ['unanchored_unlinked_keep',   false, false, false],
        ]

        TEST_CUTADAPT_GRID(
            primers_found_reads_ch,
            F_primer_ch,
            R_primer_ch,
            params.min_cutadapt_len,
            channel.fromList(cutadapt_combos)
        )

        cutadapt_grid_ch = TEST_CUTADAPT_GRID.out.grid_row
            .collect()
            .ifEmpty(file('NO_FILE_CUTADAPT')) // Emit a placeholder file if no samples had primers detected

        // Pick the best discard and keep combos from the cutadapt grid stats
        TEST_PICK_BEST_CUTADAPT(
            TEST_CUTADAPT_GRID.out.grid_row.collect()
        )
        best_cutadapt_combos_ch = TEST_PICK_BEST_CUTADAPT.out.summary
            .ifEmpty(file('NO_FILE_BEST_CUTADAPT')) // Emit a placeholder file if no cutadapt grid results

        // Read the best combo labels as value channels for filtering trimmed reads
        best_discard_ch = TEST_PICK_BEST_CUTADAPT.out.best_discard
            .splitText()
            .map { it.trim() }
            
        best_keep_ch = TEST_PICK_BEST_CUTADAPT.out.best_keep
            .splitText()
            .map { it.trim() }

        // Mix the strategies together into a single flat stream of valid strategies
        valid_strategies_ch = best_discard_ch.mix(best_keep_ch)

        // Filter trimmed reads to ONLY matching valid combos.
        // .combine() will evaluate every sample read against every active strategy.
        best_trimmed_reads_ch = TEST_CUTADAPT_GRID.out.trimmed_reads
            .combine(valid_strategies_ch)
            .filter { sample_id, reads, strategy ->
                reads[0].name.contains("_${strategy}_")
            } 
    }
    else {
        // If not running primer detection, emit placeholder files for summary, decision, and cutadapt to use in the report
        primer_summary_ch = channel.of(file('NO_FILE_PRIMER_SUMMARY'))
        primer_decision_ch = channel.of(file('NO_FILE_PRIMER_DECISION'))
        cutadapt_grid_ch = channel.of(file('NO_FILE_CUTADAPT'))
        best_cutadapt_combos_ch = channel.of(file('NO_FILE_BEST_CUTADAPT'))
        best_trimmed_reads_ch = channel.empty()
    }

    selected_isPaired_ch = selected_reads_ch
        .map { sample_id, reads, isPaired -> isPaired }
        .first()
    
    // ─────────────────────────────────────────────────────────────────────────
    // Step 3 — filterAndTrim grid (filter and full levels)
    // For 'filter': runs directly on selected reads, no primer check.
    // For 'full': runs regardless of primer decision — if no primers were found
    // filterAndTrim still runs on selected reads instead of trimmed reads.
    // ─────────────────────────────────────────────────────────────────────────
    if (test_level == 'filter' || test_level == 'full') {

        // Build sample manifest for TEST_FILTER_GRID from the selected sample set
        // so filterAndTrim only expects reads staged for those selected samples.
        selected_sample_ids_file_ch = selected_reads_ch
            .map { sample_id, reads, isPaired -> sample_id.toString() }
            .collectFile(name: 'selected_sample_ids.txt', newLine: true)

        samples_ch = selected_sample_ids_file_ch
            .combine(selected_isPaired_ch)
            .map { sample_IDs_file, isPaired -> tuple(sample_IDs_file, isPaired) }
            .first()

        // Define the raw fallback channel
        raw_fallback_ch = selected_reads_ch
            .map { sample_id, reads, isPaired -> reads instanceof List ? reads : [reads] }
            .collect()
            .map { reads_list -> tuple('raw', reads_list.flatten()) }

        // Define the processed trimmed reads channel 
        trimmed_reads_processed_ch = best_trimmed_reads_ch
            .map { sample_id, reads, strategy -> tuple(strategy, reads instanceof List ? reads : [reads]) }
            .groupTuple()
            .map { strategy, reads_list -> tuple(strategy, reads_list.flatten()) }

        // Read the primer decision text value safely into a value channel
        primers_found_ch = (test_level == 'full')
            ? TEST_PARSE_PRIMER_RESULTS.out.decision
                .map { file -> file.text.trim() == 'primers_found' }
            : channel.of(false)
        
        // Combine channels with the decision stream and filter them 
        trimmed_path_ch = trimmed_reads_processed_ch
            .combine(primers_found_ch)
            .filter { label, reads, primers_found -> test_level == 'full' && primers_found }
            .map { label, reads, primers_found -> tuple(label, reads) }

        fallback_path_ch = raw_fallback_ch
            .combine(primers_found_ch)
            .filter { label, reads, primers_found -> test_level != 'full' || !primers_found }
            .map { label, reads, primers_found -> tuple(label, reads) }

        // Since only one of the two channels will have data for a given sample set, mix them into a single channel of [label, reads] tuples for the filter grid
        filter_reads_ch = trimmed_path_ch.mix(fallback_path_ch)

        // Auto-pick truncLen from raw MultiQC per-base quality data
        // This runs regardless of whether the user supplied truncLen values, since
        // the auto value is always included in the grid as an additional combo.
        TEST_PICK_TRUNCLEN(
            multiqc_data_ch,
            selected_isPaired_ch
        )

        // Parse auto-picked truncLen values from the output file into a channel of [lt, rt] pairs
        auto_trunclen_ch = TEST_PICK_TRUNCLEN.out
            .map { file ->
                def (auto_lt, auto_rt) = file.text.trim().tokenize('\t').collect { it.toInteger() }
                [auto_lt, auto_rt]
            }

        auto_lt_ch = auto_trunclen_ch.map { it[0] }
        auto_rt_ch = auto_trunclen_ch.map { it[1] }

        // ── Build filter combo grid ─────────────────────────────────────────────────────────────────────────────────
        // The grid always includes:
        //   - 0/0 (no truncation, source = "default") — always included as a baseline, unless
        //         auto picked it (source = "auto") or user supplied it (source = "user") or both (source = "both")
        //   - auto-picked truncLen (source = "auto") — from PICK_TRUNCLEN above, unless
        //         user also supplied it (source = "both")
        //   - user-supplied truncLen values (source = "user") — from params, optional
        // All truncLen pairs are crossed with all maxEE combos.
        // Each combo is a 5-element list: [lt, rt, lm, rm, source].
        // User truncLen params default to null (no user values); set them in
        // nextflow.config or via --test_trunc_left / --test_trunc_right to add extras.
        filter_combos_ch = auto_trunclen_ch.map { auto_lt, auto_rt ->

            // User-supplied truncLen values — empty list if not provided
            def lt_vals = params.test_trunc_left
                ? params.test_trunc_left.tokenize(',').collect  { it.trim().toInteger() }
                : []
            def rt_vals = params.test_trunc_right
                ? params.test_trunc_right.tokenize(',').collect { it.trim().toInteger() }
                : []
            def lm_vals = (params.test_maxEE_left  ?: "1,2").tokenize(',').collect { it.trim().toFloat() }
            def rm_vals = (params.test_maxEE_right ?: "1,2").tokenize(',').collect { it.trim().toFloat() }

            // Build deduplicated set of truncLen pairs to run, preserving source info.
            // Use a map keyed by "lt_rt" to avoid running duplicate combos.
            def trunc_pairs = [:]

            // Determine whether user supplied 0/0 explicitly
            def user_supplied_zero = lt_vals.contains(0) && rt_vals.contains(0)

            // Determine whether auto picked 0/0
            def auto_picked_zero = (auto_lt == 0 && auto_rt == 0)

            // Always include 0/0 — label depends on user/auto overlap
            def zero_source = (user_supplied_zero && auto_picked_zero) ? "both" :
                            user_supplied_zero                        ? "user" :
                            auto_picked_zero                          ? "auto" : "default"
            trunc_pairs["0_0"] = [lt: 0, rt: 0, source: zero_source]

            // Add auto-picked values (skip if 0/0, already handled above)
            def auto_key = "${auto_lt}_${auto_rt}"
            if (auto_key != "0_0") {
                // Check if user also supplied this value
                def auto_source = (lt_vals.contains(auto_lt) && rt_vals.contains(auto_rt)) ? "both" : "auto"
                trunc_pairs[auto_key] = [lt: auto_lt, rt: auto_rt, source: auto_source]
            }

            // Add user-supplied values not already handled
            lt_vals.eachWithIndex { lt, i ->
                def rt  = rt_vals.size() > i ? rt_vals[i] : 0
                def key = "${lt}_${rt}"
                if (key != "0_0" && key != auto_key) {
                    // Distinct user value that doesn't overlap with auto or 0/0
                    trunc_pairs[key] = [lt: lt, rt: rt, source: "user"]
                }
            }

            // Build full Cartesian product of truncLen pairs × maxEE combos
            def combos = []
            trunc_pairs.values().each { pair ->
                lm_vals.each { lm ->
                    rm_vals.each { rm ->
                        combos << [pair.lt, pair.rt, lm, rm, pair.source]
                    }
                }
            }
            return combos
        }

        filter_jobs_ch = filter_reads_ch
            .combine(filter_combos_ch.flatMap { it }) 
            .map { label, reads, lt, rt, lm, rm, source ->
                tuple(label, reads, [lt, rt, lm, rm, source])
            }

        // raw_read_counts.tsv from TEST_PICK_BEST_CUTADAPT if cutadapt ran, else a placeholder.
        // .first() makes it a value channel so it broadcasts to every filter-grid job.
        raw_read_counts_ch = (test_level == 'filter'
            ? channel.value(file('NO_FILE'))
            : TEST_PICK_BEST_CUTADAPT.out.raw_read_counts
                .ifEmpty(file('NO_FILE'))
                .first()
        )

        TEST_FILTER_GRID(
            samples_ch,
            filter_jobs_ch,
            raw_read_counts_ch
        )
        filter_grid_ch = TEST_FILTER_GRID.out.read_counts_tracking.collect()
    }
    else {
        auto_lt_ch = channel.of('none')
        auto_rt_ch = channel.of('none')
        filter_grid_ch = channel.of(file('NO_FILE_FILTER'))
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Final report
    // ─────────────────────────────────────────────────────────────────────────
    acc_label = params.accession ?: params.input_file ?: 'unknown'

    TEST_REPORT(
        sample_quality_rank_ch,
        primer_summary_ch,
        primer_decision_ch,
        cutadapt_grid_ch,
        best_cutadapt_combos_ch,
        filter_grid_ch,
        test_level,
        F_primer_ch,
        R_primer_ch,
        params.target_region,
        acc_label,
        multiqc_data_ch,
        selected_isPaired_ch,
        auto_lt_ch,
        auto_rt_ch
    )

    emit:
    report = TEST_REPORT.out.report
}