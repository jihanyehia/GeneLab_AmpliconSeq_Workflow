#!/usr/bin/env Rscript
#
# Run DADA2::filterAndTrim for a single truncLen / maxEE combination,
# where truncLen values come from:
#   - auto-picked from raw MultiQC per-base quality data (Q20 threshold)
#   - user-supplied values
#   - 0/0 (no truncation, always included)
# and write one per-sample read-count tracking TSV for this combo.
#
# Usage:
#   Rscript test_filter_grid.R <assay_suffix> <unique-sample-IDs-file> <is_paired> <cutadapt_label>
#                               <left_trunc> <right_trunc> <left_maxEE> <right_maxEE> <truncLen_source>
#                               <concatenate_reads_only> <outfile>
#
# Input reads are expected in the current working directory (Nextflow stages
# them there via the flat 'path reads' input declaration).

suppressPackageStartupMessages(library(dada2))
 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11) {
    stop(paste(
        "Usage: test_filter_grid.R <assay_suffix> <unique-sample-IDs-file> <is_paired> <cutadapt_label>",
        "<left_trunc> <right_trunc> <left_maxEE> <right_maxEE> <truncLen_source> <concatenate_reads_only> <outfile>"
    ))
}
 
assay_suffix           <- args[1]
sample_IDs_file        <- args[2]
is_paired              <- as.logical(args[3])
cutadapt_label         <- args[4]
lt                     <- as.integer(args[5])
rt                     <- as.integer(args[6])
lm                     <- as.numeric(args[7])
rm                     <- as.numeric(args[8])
truncLen_source        <- args[9]
concatenate_reads_only <- as.logical(args[10])
outfile                <- args[11]

if (is.na(lt) || is.na(rt) || is.na(lm) || is.na(rm)) {
    stop("left_trunc, right_trunc, left_maxEE, right_maxEE must all be numeric.")
}
if (is.na(is_paired)) {
    stop("is_paired must be TRUE or FALSE.")
}

if (is.na(concatenate_reads_only)) {
    stop("concatenate_reads_only must be TRUE or FALSE.")
}

reads_dir <- "."
outdir <- "filtered_reads"
dir.create(outdir, showWarnings = FALSE)


# ── Sample name reader ─────────────────────────────────────────────────────────────────────────────────────
read_sample_names <- function(sample_file) {
    # Try to read as CSV first
    tryCatch({
        first_val <- scan(sample_file, what = "character", nlines = 1, quiet = TRUE)[1]
        has_header <- first_val %in% c("sample_id", "Sample Name")
        data <- read.csv(sample_file, stringsAsFactors = FALSE, header = has_header)
        
        # Check for various possible sample ID column names
        sample_col <- NULL
        possible_cols <- c("sample_id", "Sample Name")
        
        for (col in possible_cols) {
            if (col %in% colnames(data)) {
                sample_col <- col
                break
            }
        }
        
        if (!is.null(sample_col)) {
            # Found a sample ID column - use it in original order
            return(data[[sample_col]])
        } else {
            # CSV but no recognizable sample column - try first column
            if (ncol(data) >= 1) {
                return(data[,1])
            } else {
                # Fall back to single column reading
                return(scan(sample_file, what="character"))
            }
        }
    }, error = function(e) {
        # Not a valid CSV - treat as single column file
        return(scan(sample_file, what="character"))
    })
}

sample.names <- read_sample_names(sample_IDs_file)


# ── Build file paths from sample names + suffixes (raw/trimmed) ────────────────────────────────────────────

# Detect whether reads are trimmed based on files present in working dir
is_trimmed <- length(list.files(".", pattern = "_trimmed\\.fastq\\.gz$")) > 0
if (is_trimmed) {
    input_file_R1_suffix <- paste0("_", cutadapt_label, "_R1_trimmed.fastq.gz")
    input_file_R2_suffix <- paste0("_", cutadapt_label, "_R2_trimmed.fastq.gz")
    input_file_SE_suffix   <- paste0("_", cutadapt_label, "_trimmed.fastq.gz")
} else {
    input_file_R1_suffix <- paste0(assay_suffix, "_R1_raw.fastq.gz")
    input_file_R2_suffix <- paste0(assay_suffix, "_R2_raw.fastq.gz")
    input_file_SE_suffix   <- paste0(assay_suffix, "_raw.fastq.gz")
}

if (is_paired) {
    forward_reads <- file.path(reads_dir, paste0(sample.names, input_file_R1_suffix))
    forward_filtered_reads <- file.path(outdir, paste0(sample.names, "_R1_filtered.fastq.gz"))
    names(forward_filtered_reads) <- sample.names

    reverse_reads <- file.path(reads_dir, paste0(sample.names, input_file_R2_suffix))
    reverse_filtered_reads <- file.path(outdir, paste0(sample.names, "_R2_filtered.fastq.gz"))
    names(reverse_filtered_reads) <- sample.names
} else {
    forward_reads <- file.path(reads_dir, paste0(sample.names, input_file_SE_suffix))
    forward_filtered_reads <- file.path(outdir, paste0(sample.names, "_filtered.fastq.gz"))
    names(forward_filtered_reads) <- sample.names
}


# ── filterAndTrim ──────────────────────────────────────────────────────────────────────────────────────────
filtered_out <- tryCatch({
    if (is_paired) {
        filterAndTrim(
            forward_reads, forward_filtered_reads,
            reverse_reads, reverse_filtered_reads,
            truncLen    = c(lt, rt),
            maxEE       = c(lm, rm),
            maxN        = 0,
            truncQ      = 2,
            rm.phix     = TRUE,
            compress    = TRUE,
            multithread = TRUE
        )
    } else {
        filterAndTrim(
            forward_reads, forward_filtered_reads,
            truncLen    = lt,
            maxEE       = lm,
            maxN        = 0,
            truncQ      = 2,
            rm.phix     = TRUE,
            compress    = TRUE,
            multithread = TRUE
        )
    }
}, error = function(e) {
    message("filterAndTrim error: ", conditionMessage(e))
    NULL
})

if (is.null(filtered_out)) {
    stop("filterAndTrim failed; cannot proceed with denoising.")
}

# Build per-sample input counts for reporting.
# In trimmed mode, reads.in from filterAndTrim are cutadapt-trimmed counts, so
# we recover raw input counts from the staged raw_read_counts.tsv file
# (emitted by TEST_PICK_BEST_CUTADAPT from the cutadapt grid TSVs).
trimmed_reads_by_sample <- as.numeric(filtered_out[, "reads.in"])
raw_reads_by_sample <- trimmed_reads_by_sample

if (is_trimmed) {
    raw_counts_file <- "raw_read_counts.tsv"
    if (!file.exists(raw_counts_file)) {
        stop("Trimmed reads detected but raw_read_counts.tsv was not staged.")
    }

    raw_counts <- read.table(raw_counts_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (!all(c("sample", "reads_in") %in% colnames(raw_counts))) {
        stop("raw_read_counts.tsv is missing required columns: sample, reads_in")
    }

    raw_reads_by_sample <- as.numeric(raw_counts$reads_in[match(sample.names, raw_counts$sample)])
    if (any(is.na(raw_reads_by_sample))) {
        missing_samples <- sample.names[is.na(raw_reads_by_sample)]
        stop(paste0(
            "Missing raw read counts for sample(s) in raw_read_counts.tsv: ",
            paste(missing_samples, collapse = ", ")
        ))
    }
}

# Check if any reads survived filtering, and if not, write a zero-count tracking table and exit gracefully
surviving_samples <- sample.names[file.exists(forward_filtered_reads)]
if (length(surviving_samples) == 0) {
    message("No reads passed filtering for truncLen (", lt, ", ", rt, 
            ") maxEE (", lm, ", ", rm, ") — writing zero-count tracking table.")

    if (is_trimmed) {
        count_summary_tab <- data.frame(
            sample                    = sample.names,
            cutadapt_label            = cutadapt_label,
            truncLen_source           = truncLen_source,
            left_trunc                = lt,
            right_trunc               = rt,
            left_maxEE                = lm,
            right_maxEE               = rm,
            raw_reads                 = raw_reads_by_sample,
            cutadapt_trimmed          = trimmed_reads_by_sample,
            dada2_filtered            = 0,
            dada2_denoised_F          = 0,
            dada2_denoised_R          = if (is_paired) 0 else NULL,
            dada2_merged              = if (is_paired) 0 else NULL,
            dada2_chimera_removed     = 0,
            final_perc_reads_retained = 0,
            notes                     = "All reads removed by filterAndTrim"
        )
    } else {
        count_summary_tab <- data.frame(
            sample                    = sample.names,
            cutadapt_label            = cutadapt_label,
            truncLen_source           = truncLen_source,
            left_trunc                = lt,
            right_trunc               = rt,
            left_maxEE                = lm,
            right_maxEE               = rm,
            raw_reads                 = raw_reads_by_sample,
            dada2_filtered            = 0,
            dada2_denoised_F          = 0,
            dada2_denoised_R          = if (is_paired) 0 else NULL,
            dada2_merged              = if (is_paired) 0 else NULL,
            dada2_chimera_removed     = 0,
            final_perc_reads_retained = 0,
            notes                     = "All reads removed by filterAndTrim"
        )
    }
    
    if (!is_paired) {
        count_summary_tab$dada2_denoised_R <- NULL
        count_summary_tab$dada2_merged     <- NULL
    }
    
    write.table(count_summary_tab, file = outfile, sep = "\t", 
                quote = FALSE, row.names = FALSE)
    message("Written: ", outfile)
    quit(status = 0)
}


# ── Learn errors ───────────────────────────────────────────────────────────────────────────────────────────
# Only pass filtered reads that actually exist because aggressive truncation/maxEE combos may drop all reads
# Samples with 0 output reads will have empty files that learnErrors cannot handle

fwd_exist <- forward_filtered_reads[file.exists(forward_filtered_reads)]
forward_errors <- learnErrors(fwd_exist, multithread = TRUE)
forward_seqs   <- dada(fwd_exist, err = forward_errors, pool = "pseudo", multithread = TRUE)

if (is_paired) {
    rev_exist <- reverse_filtered_reads[file.exists(reverse_filtered_reads)]
    reverse_errors <- learnErrors(rev_exist, multithread = TRUE)
    reverse_seqs   <- dada(rev_exist, err = reverse_errors, pool = "pseudo", multithread = TRUE)
}


# ── Merge (PE) or proceed directly (SE) ────────────────────────────────────────────────────────────────────
if (is_paired) {
    if (concatenate_reads_only) {
        merged_contigs <- mergePairs(
            forward_seqs, fwd_exist,
            reverse_seqs, rev_exist,
            verbose          = TRUE,
            justConcatenate  = TRUE
        )
    } else {
        merged_contigs <- mergePairs(
            forward_seqs, fwd_exist,
            reverse_seqs, rev_exist,
            verbose = TRUE
        )
    }
    seqtab <- makeSequenceTable(merged_contigs)
} else {
    seqtab <- makeSequenceTable(forward_seqs)
}


# ── Remove chimeras ────────────────────────────────────────────────────────────────────────────────────────
seqtab.nochim <- removeBimeraDenovo(
    seqtab, method = "consensus", multithread = TRUE, verbose = TRUE
)


# ── Per-sample read-count tracking table ───────────────────────────────────────────────────────────────────

# Helper function to count unique sequences in a dada2 object
getN <- function(x) sum(getUniques(x))

# Ensure consistent sample ordering for all downstream steps by using only samples that made it through DADA2
samples_in_seqtab <- rownames(seqtab.nochim)

# filtered_out rows are named by input path; align to sample names
filtered_count_tab <- data.frame(
    sample           = sample.names,
    raw_reads        = raw_reads_by_sample,
    cutadapt_trimmed = trimmed_reads_by_sample,
    reads_out        = filtered_out[, "reads.out"]
)

# Re-index filtered_count_tab against all original samples, not just survivors
filtered_count_tab_ordered <- filtered_count_tab[
    match(samples_in_seqtab, filtered_count_tab$sample), 
]

if (is_paired) {
    if (is_trimmed) {
        count_summary_tab <- data.frame(
            sample                = samples_in_seqtab,
            cutadapt_label        = cutadapt_label,
            truncLen_source       = truncLen_source,
            left_trunc            = lt,
            right_trunc           = rt,
            left_maxEE            = lm,
            right_maxEE           = rm,
            raw_reads             = filtered_count_tab_ordered$raw_reads,
            cutadapt_trimmed      = filtered_count_tab_ordered$cutadapt_trimmed,
            dada2_filtered        = filtered_count_tab_ordered$reads_out,
            dada2_denoised_F      = sapply(forward_seqs, getN),
            dada2_denoised_R      = sapply(reverse_seqs, getN),
            dada2_merged          = rowSums(seqtab)[samples_in_seqtab],
            dada2_chimera_removed = rowSums(seqtab.nochim),
            final_perc_reads_retained = ifelse(
                filtered_count_tab_ordered$raw_reads > 0,
                round(rowSums(seqtab.nochim) / filtered_count_tab_ordered$raw_reads * 100, 1),
                0
            )
        )
    } else {
        count_summary_tab <- data.frame(
            sample                = samples_in_seqtab,
            cutadapt_label        = cutadapt_label,
            truncLen_source       = truncLen_source,
            left_trunc            = lt,
            right_trunc           = rt,
            left_maxEE            = lm,
            right_maxEE           = rm,
            raw_reads             = filtered_count_tab_ordered$raw_reads,
            dada2_filtered        = filtered_count_tab_ordered$reads_out,
            dada2_denoised_F      = sapply(forward_seqs, getN),
            dada2_denoised_R      = sapply(reverse_seqs, getN),
            dada2_merged          = rowSums(seqtab)[samples_in_seqtab],
            dada2_chimera_removed = rowSums(seqtab.nochim),
            final_perc_reads_retained = ifelse(
                filtered_count_tab_ordered$raw_reads > 0,
                round(rowSums(seqtab.nochim) / filtered_count_tab_ordered$raw_reads * 100, 1),
                0
            )
        )
    }
} else {
    if (is_trimmed) {
        count_summary_tab <- data.frame(
            sample                = samples_in_seqtab,
            cutadapt_label        = cutadapt_label,
            truncLen_source       = truncLen_source,
            left_trunc            = lt,
            right_trunc           = rt,
            left_maxEE            = lm,
            right_maxEE           = rm,
            raw_reads             = filtered_count_tab_ordered$raw_reads,
            cutadapt_trimmed      = filtered_count_tab_ordered$cutadapt_trimmed,
            dada2_filtered        = filtered_count_tab_ordered$reads_out,
            dada2_denoised        = sapply(forward_seqs, getN),
            dada2_chimera_removed = rowSums(seqtab.nochim),
            final_perc_reads_retained = ifelse(
                filtered_count_tab_ordered$raw_reads > 0,
                round(rowSums(seqtab.nochim) / filtered_count_tab_ordered$raw_reads * 100, 1),
                0
            )
        )
    } else {
        count_summary_tab <- data.frame(
            sample                = samples_in_seqtab,
            cutadapt_label        = cutadapt_label,
            truncLen_source       = truncLen_source,
            left_trunc            = lt,
            right_trunc           = rt,
            left_maxEE            = lm,
            right_maxEE           = rm,
            raw_reads             = filtered_count_tab_ordered$raw_reads,
            dada2_filtered        = filtered_count_tab_ordered$reads_out,
            dada2_denoised        = sapply(forward_seqs, getN),
            dada2_chimera_removed = rowSums(seqtab.nochim),
            final_perc_reads_retained = ifelse(
                filtered_count_tab_ordered$raw_reads > 0,
                round(rowSums(seqtab.nochim) / filtered_count_tab_ordered$raw_reads * 100, 1),
                0
            )
        )
    }
}

write.table(count_summary_tab, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
message("Written: ", outfile)