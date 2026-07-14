#!/usr/bin/env python3
"""
Summarise primer detection results and decide whether primers are present.

Reads all *_primer_hits.tsv files in the current directory (staged by Nextflow),
computes the mean of the per-sample maximum primer hit rate across all four
orientations (F, F_rc, R, R_rc) on R1 reads, and writes a decision file
that the subworkflow uses to branch into cutadapt grid vs direct filterAndTrim.

Decision logic:
    For each sample, take the highest hit rate across all four primer orientations
    on R1. Average this across all samples. If the mean exceeds the threshold,
    primers are called present.

    This approach is orientation-agnostic — it correctly handles cases where
    primers appear as readthrough RC hits rather than 5' forward hits, which
    is common for short amplicons where the read length exceeds the insert size.

Output files:
    primer_detection_summary.tsv  — all primer hit rows consolidated
    primer_decision.txt           — single line: "primers_found" or "primers_not_found"
"""

import argparse
import csv
import glob
import sys

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--threshold',       type=float, default=5.0,
                        help='Min mean max-primer-pct on R1 to call primers present (default: 5.0)')
    parser.add_argument('--output-summary',  required=True,
                        help='Path for per-sample summary TSV')
    parser.add_argument('--output-decision', required=True,
                        help='Path for decision file (primers_found / primers_not_found)')
    args = parser.parse_args()

    PRIMER_COLS = ('F_primer_pct', 'F_primer_rc_pct', 'R_primer_pct', 'R_primer_rc_pct')

    # ── load all primer hit TSVs staged into CWD ─────────────────────────────────
    all_rows = []
    for f in sorted(glob.glob('*_primer_hits.tsv')):
        try:
            with open(f) as fh:
                all_rows.extend(csv.DictReader(fh, delimiter='\t'))
        except Exception as exc:
            print(f"Warning: could not read {f}: {exc}", file=sys.stderr)

    if not all_rows:
        sys.exit("ERROR: no *_primer_hits.tsv files found in working directory.")

    # ── compute decision ──────────────────────────────────────────────────────────
    # Use R1 rows for the decision; fall back to all rows for single-end.
    r1_rows = [r for r in all_rows if r.get('read') == 'R1'] or all_rows

    max_pcts = []
    for r in r1_rows:
        try:
            vals = [float(r[col]) for col in PRIMER_COLS if r.get(col)]
            if vals:
                max_pcts.append(max(vals))
        except (KeyError, ValueError):
            pass

    if not max_pcts:
        print("Warning: could not parse primer hit values; assuming primers not found.",
            file=sys.stderr)
        mean_max_pct = 0.0
        decision   = "primers_not_found"
    else:
        mean_max_pct = sum(max_pcts) / len(max_pcts)
        decision   = "primers_found" if mean_max_pct >= args.threshold else "primers_not_found"

    print(f"Mean max primer pct on R1: {mean_max_pct:.1f}% → {decision}")

    # ── write outputs ─────────────────────────────────────────────────────────────
    summary_cols = list(all_rows[0].keys()) if all_rows else []
    with open(args.output_summary, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=summary_cols, delimiter='\t',
                        extrasaction='ignore')
        w.writeheader()
        w.writerows(all_rows)

    with open(args.output_decision, 'w') as fh:
        fh.write(decision + '\n')

    print(f"Written: {args.output_summary}, {args.output_decision}")

if __name__ == "__main__":
    main()