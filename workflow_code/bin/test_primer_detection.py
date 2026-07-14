#!/usr/bin/env python
"""
Scan the first N reads of a FASTQ(.gz) file for forward and reverse primers
and their reverse complements. Writes a per-sample TSV of hit counts and rates.

Output file: <sample_id>_primer_hits.tsv
Columns: sample, read, total_checked,
         F_primer_hits, F_primer_pct, F_primer_rc_hits, F_primer_rc_pct,
         R_primer_hits, R_primer_pct, R_primer_rc_hits, R_primer_rc_pct
"""

import argparse
import re
import gzip
import sys
import csv

def rc(seq):
    comp = str.maketrans('ACGTNRYWSKMBDHV', 'TGCANYRWSMKVHDB')
    return seq.translate(comp)[::-1]

def primer_regex(primer):
    """Convert a primer sequence with IUPAC codes to a regex pattern."""
    IUPAC = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }

    return re.compile(''.join(IUPAC[b] for b in primer.upper()))

def count_hits(filepath, primers, n):
    """Return (hits_dict, total_reads_checked)."""
    openfn = gzip.open if filepath.endswith('.gz') else open
    hits   = {p: 0 for p in primers}
    total  = 0
    try:
        with openfn(filepath, 'rt') as fh:
            while total < n:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().strip().upper()
                fh.readline()   # +
                fh.readline()   # qual
                total += 1
                for name, pattern in primers.items():
                    if pattern.search(seq):
                        hits[name] += 1
    except Exception as e:
        print(f"Warning: could not read {filepath}: {e}", file=sys.stderr)
    return hits, total

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-id', required=True, help='Sample identifier for reporting')
    parser.add_argument('--r1', required=True, help='Path to R1 FASTQ')
    parser.add_argument('--r2', help='Path to R2 FASTQ (optional)')
    parser.add_argument('--f-primer', required=True, help='Forward primer sequence')
    parser.add_argument('--r-primer', required=True, help='Reverse primer sequence')
    parser.add_argument('--n-reads', type=int, default=1000, help='Number of reads to check from each FASTQ')
    args = parser.parse_args()

    f_primer = args.f_primer.upper()
    r_primer = args.r_primer.upper()

    primers_to_check = {
        'F_primer':    primer_regex(f_primer),
        'F_primer_rc': primer_regex(rc(f_primer)),
        'R_primer':    primer_regex(r_primer),
        'R_primer_rc': primer_regex(rc(r_primer)),
    }
    hits_r1, total_r1 = count_hits(args.r1, primers_to_check, args.n_reads)

    hits_r2, total_r2 = {}, 0
    if args.r2:
        hits_r2, total_r2 = count_hits(args.r2, primers_to_check, args.n_reads)

    outfile = f"{args.sample_id}_primer_hits.tsv"
    with open(outfile, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['sample','read','total_checked',
                    'F_primer_hits','F_primer_pct',
                    'F_primer_rc_hits','F_primer_rc_pct',
                    'R_primer_hits','R_primer_pct',
                    'R_primer_rc_hits','R_primer_rc_pct'])
        
        def write_row(label, hits, total) -> None:
            
            def fmt(key):
                h = hits.get(key,0)
                pct = round(100*h/total,1) if total else 0
                return h, pct
            
            fh, fp   = fmt('F_primer')
            frh, frp = fmt('F_primer_rc')
            rh, rp   = fmt('R_primer')
            rrh, rrp = fmt('R_primer_rc')
            w.writerow([args.sample_id, label, total,
                        fh, fp, frh, frp, rh, rp, rrh, rrp])
        write_row('R1', hits_r1, total_r1)
        if args.r2:
            write_row('R2', hits_r2, total_r2)
    
    print(f"Written: {outfile}")

if __name__ == "__main__":
    main()