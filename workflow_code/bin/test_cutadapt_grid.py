#!/usr/bin/env python
"""
Parse a cutadapt --json output file and write one TSV row summarising the run.
Called after cutadapt finishes in the TEST_CUTADAPT_GRID process.

Output file: <sample_id>_<combo_label>_cutadapt_grid.tsv
Columns: sample, combo, anchored, linked, discard_untrimmed,
         reads_in, reads_out, pct_retained, pct_with_adapter
"""

import argparse
import csv
import json

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--sample-id", required=True, help='Sample identifier for reporting')
    parser.add_argument("--combo-label", required=True, help='Label for the primer combination')
    parser.add_argument("--anchored", required=True, help='Anchored primer status')
    parser.add_argument("--linked", required=True, help='Linked primer status')
    parser.add_argument("--discard", required=True, help='Discard untrimmed reads status')
    parser.add_argument("--json-file", required=True, help='Path to the cutadapt JSON output file')

    args = parser.parse_args()

    try:

        with open(args.json_file) as fh:
            data = json.load(fh)

        read_counts = data.get("read_counts", {})
        reads_in = read_counts.get("input", 0)
        reads_out = read_counts.get("output", 0)
        pct_retained = round(100 * reads_out / reads_in, 1) if reads_in else 0

        read1_with_adapter = read_counts.get("read1_with_adapter") or 0
        read2_with_adapter = read_counts.get("read2_with_adapter") # May be None for single-end data
        pct_r1_with_adapter = round(100 * read1_with_adapter / reads_in, 1) if reads_in else 0
        is_paired = read2_with_adapter is not None
        pct_r2_with_adapter = round(100 * read2_with_adapter / reads_in, 1) if reads_in and is_paired else 'N/A'

    except Exception as exc:
        print(f"Warning: could not parse {args.json_file}: {exc}", file=sys.stderr)
        reads_in = reads_out = with_adapter = 'ERROR'
        pct_retained = pct_with_adapter = 'ERROR'

    outfile = f"{args.sample_id}_{args.combo_label}_cutadapt_grid.tsv"
    with open(outfile, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow([
        'sample', 'combo', 'anchored', 'linked', 'discard_untrimmed',
        'reads_in', 'reads_out', 'pct_retained', 'pct_r1_with_adapter', 'pct_r2_with_adapter'
        ])
        w.writerow([
            args.sample_id, args.combo_label, args.anchored, args.linked, args.discard,
            reads_in, reads_out, pct_retained, pct_r1_with_adapter, pct_r2_with_adapter
        ])
    print(f"Written: {outfile}")


if __name__ == "__main__":
    main()