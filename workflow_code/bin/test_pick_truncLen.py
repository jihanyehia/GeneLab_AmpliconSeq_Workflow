#!/usr/bin/env python3
"""
Pick reasonable DADA2 truncLen values from raw MultiQC per-base sequence
quality data using a Q20 drop threshold.

Usage:
    pick_trunclen.py --multiqc-json <path> --is-paired <true|false>
                     --output <path>

Output file contains two tab-separated integers: auto_lt  auto_rt
For SE data, auto_rt is always 0.
"""

import argparse
import json


def pick_trunclen(lines, direction, min_qual=20):
    """
    Given a list of MultiQC per-base quality lines, pick a truncLen for the
    specified direction (R1 or R2) based on where average quality drops below
    min_qual.
    """
    dir_lines = [l for l in lines if l['name'].endswith('_' + direction)]

    if not dir_lines:
        print(f"No {direction} lines found in MultiQC JSON; returning 0.")
        return 0

    # Collect all positions across all samples
    all_positions = sorted(set(p[0] for l in dir_lines for p in l['pairs']))

    # Compute median quality per position across all samples
    median_quals = []
    for pos in all_positions:
        vals = sorted(p[1] for l in dir_lines for p in l['pairs'] if p[0] == pos)
        if vals:
            n = len(vals)
            median = (vals[n // 2] if n % 2 == 1
                      else (vals[n // 2 - 1] + vals[n // 2]) / 2)
            median_quals.append(median)
        else:
            median_quals.append(None)

    # Find first position where median quality drops below min_qual
    drop_idx = next(
        (i for i, q in enumerate(median_quals) if q is not None and q < min_qual),
        None
    )

    if drop_idx is None:
        trunc = 0
        print(f"{direction}: quality never drops below Q{min_qual}; "
              f"returning 0 (no truncation recommended).")
    elif drop_idx == 0:
        trunc = 0
        print(f"{direction}: quality drops below Q{min_qual} at first position; "
              f"returning 0 (truncation doesn't help).")
    else:
        trunc = int(all_positions[drop_idx - 1])
        print(f"{direction}: quality drops below Q{min_qual} at position "
              f"{all_positions[drop_idx]}; truncLen set to {trunc}.")

    # Round down to nearest 10, minimum 50
    trunc = 0 if trunc == 0 else max(50, (trunc // 10) * 10)
    print(f"{direction}: final auto truncLen (rounded) = {trunc}.")
    return trunc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--multiqc-json", required=True,
                        help="Path to multiqc_data.json")
    parser.add_argument("--is-paired",    required=True,
                        help="Whether data is paired-end (true/false)")
    parser.add_argument("--output",       required=True,
                        help="Output file path for auto truncLen values")
    args = parser.parse_args()

    is_paired = args.is_paired.strip().lower() == "true"

    with open(args.multiqc_json) as f:
        d = json.load(f)

    lines = (d["report_plot_data"]
              ["fastqc_per_base_sequence_quality_plot"]
              ["datasets"][0]
              ["lines"])

    auto_lt = pick_trunclen(lines, "R1", min_qual=20)
    auto_rt = pick_trunclen(lines, "R2", min_qual=20) if is_paired else 0

    with open(args.output, "w") as f:
        f.write(f"{auto_lt}\t{auto_rt}\n")

    print(f"Written: {args.output} (auto_lt={auto_lt}, auto_rt={auto_rt})")


if __name__ == "__main__":
    main()