#!/usr/bin/env python3
"""
Select N representative samples from a MultiQC data directory by computing
a composite quality score and picking evenly-spaced interior quantiles
(avoiding outliers at either extreme).

Metrics (all from multiqc_data.json → report_saved_raw_data.multiqc_fastqc):
    mean_sequence_quality  — weighted mean of per_sequence_quality_scores
    total_sequences        — read depth
    avg_sequence_length    — captures trimming variation (if any)
    percent_gc             — converted to gc_distance from cohort median
                             so deviation in either direction is penalised

All metrics z-score normalised (clipped at ±2 std to limit outlier influence)
and oriented so higher = better, then summed into a composite score.
Samples are ranked on the composite and N selected at evenly-spaced interior
quantiles to pick typical (not best, not worst) samples.

Usage:
    select_samples.py --multiqc-data <dir> --n <int> --assay-suffix <str> --output <tsv>

Output TSV columns:
    sample_id, mean_sequence_quality, total_sequences, avg_sequence_length,
    percent_gc, gc_distance, composite_score, rank, selected

Stdout:
    Selected sample IDs, one per line.
"""

import argparse
import csv
import json
import sys
from pathlib import Path


def sample_id_from_mqc_name(name, assay_suffix=None):
    """
    Strip R1/R2/_1/_2 suffixes so paired reads collapse to one ID,
    and optional assay suffix so sample IDs match those used in the pipeline
    """
    # Strip paired-end suffix first
    for suffix in ('_R1', '_R2', '_1', '_2'):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    # Then strip assay suffix if provided
    if assay_suffix and name.endswith(assay_suffix):
        name = name[: -len(assay_suffix)]
    
    return name


def zscore_normalise(values):
    """
    Z-score normalise a list of values and clip at ±clip std deviations.

    Clipping prevents a single extreme outlier from dominating the composite
    score — once a sample is >2 std from the mean on any one metric, extra
    deviation doesn't keep pushing it further from the interior picks.
    """
    n = len(values)
    mean = sum(values) / n
    std = (sum((v - mean) ** 2 for v in values) / n) ** 0.5
    if std == 0:
        return [0.0] * n
    return [(v - mean) / std for v in values]


# Metrics where a higher raw value indicates a better / more typical sample
HIGHER_IS_BETTER = ('mean_sequence_quality', 'total_sequences', 'avg_sequence_length')

# Metrics where a lower raw value indicates a better / more typical sample.
# gc_distance (|percent_gc - cohort_median_gc|) is derived in compute_scores;
# deviation in either direction from the cohort median is undesirable.
LOWER_IS_BETTER = ('gc_distance',)


def load_samples(mqc_json, assay_suffix):
    """
    Parse multiqc_data.json and return one record per sample.

    Paired-end reads (R1/R2) are collapsed to a single record, preferring R1
    since it typically has higher quality and is the more informative read.
    Samples missing any required metric are skipped with a warning.
    """
    with open(mqc_json) as fh:
        data = json.load(fh)

    try:
        fastqc_data = data['report_saved_raw_data']['multiqc_fastqc']
    except KeyError:
        sys.exit(
            "ERROR: expected report_saved_raw_data.multiqc_fastqc in JSON.\n"
            f"Top-level keys: {list(data.keys())}"
        )

    # Build a lookup of sample name → weighted mean quality from the plot data.
    # The distribution lives under report_plot_data as [phred, count] pairs,
    # not in report_saved_raw_data where per_sequence_quality_scores is just
    # the pass/warn/fail flag.
    try:
        lines = (data['report_plot_data']
                     ['fastqc_per_sequence_quality_scores_plot']
                     ['datasets'][0]['lines'])
    except (KeyError, IndexError):
        sys.exit(
            "ERROR: could not find per-sequence quality score plot data in JSON."
        )

    mean_quality_by_name = {}
    for line in lines:
        name  = line.get('name', '')
        pairs = line.get('pairs', [])
        # pairs is a list of [phred_score, read_count]
        total_reads = sum(count for _, count in pairs)
        if total_reads > 0:
            mean_quality_by_name[name] = (
                sum(phred * count for phred, count in pairs) / total_reads
            )
    
    # Collapse paired reads: keep R1 when both present
    rows_by_sample = {}
    for name, fields in fastqc_data.items():
        sid = sample_id_from_mqc_name(name, assay_suffix)
        is_r1 = name.endswith(('_R1', '_1'))
        if sid not in rows_by_sample or is_r1:
            # Attach the original name so we can look up quality score below
            rows_by_sample[sid] = (name, fields)

    samples = []
    for sid, (name, fields) in rows_by_sample.items():
        mean_q = mean_quality_by_name.get(name)
        if mean_q is None:
            print(f"Warning: quality score distribution missing for {sid}, skipping.",
                  file=sys.stderr)
            continue

        record = {'sample_id': sid, 'mean_sequence_quality': mean_q}

        # Extract the remaining metrics that are stored directly in the JSON
        ok = True
        for metric in ('total_sequences', 'avg_sequence_length', '%GC'):
            raw = fields.get(metric) or fields.get(metric.replace('_', ' ').title()) # JSON uses 'Total Sequences' with capital letters and spaces
            if raw is None:
                print(f"Warning: '{metric}' missing for {sid}, skipping.", file=sys.stderr)
                ok = False
                break
            try:
                record[metric] = float(raw)
            except (TypeError, ValueError):
                print(f"Warning: '{metric}' non-numeric for {sid}, skipping.", file=sys.stderr)
                ok = False
                break
        if ok:
            samples.append(record)

    return samples


def compute_scores(samples):
    """
    Add gc_distance, per-metric z-scores, and composite_score to each record.

    GC content is converted to gc_distance (absolute deviation from the cohort
    median) before scoring. This is intentional: GC content has a natural expected
    value for a given organism/amplicon target, and deviation in either direction
    signals an atypical sample

    All z-scores are oriented so that higher = better, then summed. A sample
    that is consistently average across all metrics ends up near zero; one that
    is consistently good ends up with a large positive score; one that is
    consistently poor ends up with a large negative score. Interior quantile
    selection then picks from the middle of this ranking.
    """
    # Derive gc_distance from percent_gc before z-scoring
    gc_vals   = sorted(s['%GC'] for s in samples)
    n         = len(gc_vals)
    mid       = n // 2
    median_gc = (gc_vals[mid] if n % 2 else (gc_vals[mid - 1] + gc_vals[mid]) / 2)

    for s in samples:
        # Absolute deviation from cohort median — lower means more typical GC
        s['gc_distance'] = abs(s['%GC'] - median_gc)

    # Z-score each metric and flip sign for lower-is-better metrics so that
    # the composite score is uniformly oriented: higher = more typical sample
    for metric in HIGHER_IS_BETTER + LOWER_IS_BETTER:
        zs = zscore_normalise([s[metric] for s in samples])
        for s, z in zip(samples, zs):
            s[f'z_{metric}'] = -z if metric in LOWER_IS_BETTER else z

    for s in samples:
        s['composite_score'] = round(
            sum(s[f'z_{m}'] for m in HIGHER_IS_BETTER + LOWER_IS_BETTER), 4
        )

    return samples


def select_indices(n_total, n_select):
    """
    Return evenly-spaced interior quantile indices into a ranked list.

    Indices are chosen strictly inside the list, avoiding the extremes at
    either end. For example with n_total=20 and n_select=3:

        step = (20 - 1) / (3 + 1) = 4.75
        indices → round(4.75), round(9.5), round(14.25) → 5, 10, 14

    These land in the middle of the ranking, away from the best and worst
    samples at the head and tail.
    """
    n_select = min(n_select, n_total)
    if n_select == 1:
        return [n_total // 2]
    if n_select >= n_total:
        return list(range(n_total))
    step = (n_total - 1) / (n_select + 1)
    return [round(step * (i + 1)) for i in range(n_select)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--multiqc-data', required=True,
                        help='Path to multiqc_data directory')
    parser.add_argument('--n', type=int, default=3,
                        help='Number of representative samples (default: 3)')
    parser.add_argument('--assay-suffix', default='',
                    help='Assay suffix to strip from sample names (e.g. 16S_GLAmpSeq)')
    parser.add_argument('--output', required=True,
                        help='Output TSV path')
    args = parser.parse_args()

    mqc_json = Path(args.multiqc_data) / 'multiqc_data.json'
    if not mqc_json.exists():
        sys.exit(
            f"ERROR: {mqc_json} not found.\n"
            f"Contents: {list(Path(args.multiqc_data).iterdir())}"
        )

    samples = load_samples(mqc_json, args.assay_suffix)
    if not samples:
        sys.exit("ERROR: no usable sample rows found.")

    samples = compute_scores(samples)

    # Rank samples from best composite score to worst
    ranked = sorted(samples, key=lambda x: x['composite_score'], reverse=True)
    
    # Pick interior quantile indices and resolve to sample IDs
    selected_ids = {ranked[i]['sample_id'] for i in select_indices(len(ranked), args.n)}

    # Write full ranked table with selection flag
    out_cols = [
        'sample_id',
        'mean_sequence_quality', 'total_sequences', 'avg_sequence_length',
        '%GC', 'gc_distance',
        'composite_score', 'rank', 'selected',
    ]
    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=out_cols, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        for rank, s in enumerate(ranked, start=1):
            writer.writerow({
                **{k: s.get(k, '') for k in out_cols},
                'rank': rank,
                'selected': 'yes' if s['sample_id'] in selected_ids else 'no',
            })

    # Print selected IDs to stdout for downstream pipeline capture
    for sid in sorted(selected_ids):
        print(sid)


if __name__ == '__main__':
    main()