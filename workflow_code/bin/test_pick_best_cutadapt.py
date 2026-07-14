#!/usr/bin/env python3
"""
Pick the best cutadapt combo(s) from the grid results for downstream filtering.

Selection logic:
    Valid combos:
        valid_discard: mean_retained > 0%
        valid_keep:    mean_retained < 100% AND mean_max_adapter > 0%

    Ranking:
        Keep combos are ranked by:
            1. max(pct_r1_with_adapter, pct_r2_with_adapter) descending (prefer combos that actually trimmed something)
            2. pct_retained descending
            3. anchored > unanchored
            4. linked > unlinked

        Discard combos are ranked by:
            1. pct_retained descending
            2. max(pct_r1_with_adapter, pct_r2_with_adapter) descending
            3. anchored > unanchored
            4. linked > unlinked

    Closeness tiebreaker (CLOSENESS_THRESHOLD = 2%):
        If multiple combos are within 2% of the best on both mean_retained and
        mean_max_adapter, prefer anchored over unanchored and linked over unlinked.
        A note is added to the output when this tiebreaker is applied.

    Selection (threshold = 50%):
        best_discard >= 50% AND best_keep >= 50%:
               → best discard + best keep
        best_discard >= 50% (keep not valid or < 50%):
               → up to two best valid_discard above threshold
               → warn if no valid keep
        best_keep >= 50% (discard < 50%):
               → up to two best valid_keep above threshold
               → warn: discard had low retention
        both below 50%:
               → top two from valid_discard + valid_keep combined
               → warn: both groups had low retention

Output files:
    best_cutadapt_combos.tsv  — selected combos with mean stats and warnings
    best_discard_combo.txt    — selected discard combo label(s), one per line
    best_keep_combo.txt       — selected keep combo label(s), one per line
"""

import argparse
import csv
import glob
import sys


RETENTION_THRESHOLD = 50.0
CLOSENESS_THRESHOLD =  2.0

def parse_float(val, default=0.0):
    """Safely parse a float, returning default for NA or invalid values."""
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def anchored_score(combo_label):
    return 1 if 'anchored' in combo_label and 'unanchored' not in combo_label else 0


def linked_score(combo_label):
    return 1 if 'linked' in combo_label and 'unlinked' not in combo_label else 0


def rank_discard(entry):
    return (
        entry['mean_retained'],
        entry['mean_max_adapter'],
        anchored_score(entry['combo']),
        linked_score(entry['combo']),
    )


def rank_keep(entry):
    return (
        entry['mean_max_adapter'],
        entry['mean_retained'],
        anchored_score(entry['combo']),
        linked_score(entry['combo']),
    )


def rank_combined(entry):
    """For fallback — rank by retained then adapter regardless of group."""
    return (
        entry['mean_retained'],
        entry['mean_max_adapter'],
        anchored_score(entry['combo']),
        linked_score(entry['combo']),
    )


def pick_best_with_preference(ranked_combos):
    """
    From a ranked list, find all combos within CLOSENESS_THRESHOLD of the best
    on both mean_retained and mean_max_adapter, then prefer anchored over
    unanchored and linked over unlinked within that pool.

    Returns (winner, note) where note is None if no tiebreaker was applied.
    """
    if not ranked_combos:
        return None, None

    best          = ranked_combos[0]
    best_retained = best['mean_retained']
    best_adapter  = best['mean_max_adapter']

    # Build close pool — all combos within CLOSENESS_THRESHOLD on both metrics
    close_pool = [
        e for e in ranked_combos
        if abs(e['mean_retained']     - best_retained) <= CLOSENESS_THRESHOLD
        and abs(e['mean_max_adapter'] - best_adapter)  <= CLOSENESS_THRESHOLD
    ]

    # Within the close pool, prefer anchored then linked
    close_pool_sorted = sorted(
        close_pool,
        key=lambda e: (anchored_score(e['combo']), linked_score(e['combo'])),
        reverse=True
    )

    winner = close_pool_sorted[0]
    note   = None

    if winner['combo'] != ranked_combos[0]['combo']:
        note = (
            f"'{winner['combo']}' was selected over '{ranked_combos[0]['combo']}' "
            f"because it is within {CLOSENESS_THRESHOLD}% and is preferred "
            f"(anchored/linked over unanchored/unlinked)."
        )
        print(f"Note: {note}", file=sys.stderr)

    return winner, note


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-summary',      required=True,
                        help='Path for best combos summary TSV')
    parser.add_argument('--output-best-discard', required=True,
                        help='Path for best discard combo label file')
    parser.add_argument('--output-best-keep',    required=True,
                        help='Path for best keep combo label file')
    args = parser.parse_args()

    # ── load all cutadapt grid TSVs staged into CWD ───────────────────────────
    all_rows = []
    for f in sorted(glob.glob('*_cutadapt_grid.tsv')):
        try:
            with open(f) as fh:
                all_rows.extend(csv.DictReader(fh, delimiter='\t'))
        except Exception as exc:
            print(f"Warning: could not read {f}: {exc}", file=sys.stderr)

    if not all_rows:
        sys.exit("ERROR: no *_cutadapt_grid.tsv files found in working directory.")

    # ── aggregate stats per combo ─────────────────────────────────────────────
    # Build dict: combo_label → {discard, retained_vals, adapter_vals}
    combo_stats = {}
    for row in all_rows:
        combo   = row.get('combo', '')
        discard = row.get('discard_untrimmed', '').lower() == 'true'

        pct_retained  = parse_float(row.get('pct_retained'))
        pct_r1        = parse_float(row.get('pct_r1_with_adapter'))
        pct_r2        = parse_float(row.get('pct_r2_with_adapter'))
        max_adapter   = max(pct_r1, pct_r2)

        if combo not in combo_stats:
            combo_stats[combo] = {
                'discard':       discard,
                'anchored':      row.get('anchored', ''),
                'linked':        row.get('linked', ''),
                'retained_vals': [],
                'adapter_vals':  [],
            }
        combo_stats[combo]['retained_vals'].append(pct_retained)
        combo_stats[combo]['adapter_vals'].append(max_adapter)

    # ── compute means and split into groups ───────────────────────────────────
    discard_combos = []
    keep_combos    = []

    for combo, stats in combo_stats.items():
        n             = len(stats['retained_vals'])
        mean_retained = sum(stats['retained_vals']) / n
        mean_max_adapter  = sum(stats['adapter_vals'])  / n

        entry = {
            'combo':          combo,
            'mean_retained':  round(mean_retained, 1),
            'mean_max_adapter':   round(mean_max_adapter,  1),
        }

        if stats['discard']:
            discard_combos.append(entry)
        else:
            keep_combos.append(entry)

    # ── filter valid combos ───────────────────────────────────────────────────
    # valid_discard: retained something (> 0%)
    # valid_keep:    did real trimming (< 100% retained AND > 0% adapter hits)
    valid_discard = [combo for combo in discard_combos if combo['mean_retained'] > 0]
    valid_keep    = [combo for combo in keep_combos
                     if combo['mean_retained'] < 100.0 and combo['mean_max_adapter'] > 0]

    # ── rank each group ───────────────────────────────────────────────────────
    valid_discard_ranked = sorted(valid_discard, key=rank_discard, reverse=True)
    valid_keep_ranked    = sorted(valid_keep,    key=rank_keep,    reverse=True)

    # ── pick best with closeness tiebreaker ───────────────────────────────────
    best_discard, discard_note = (
        pick_best_with_preference(valid_discard_ranked)
        if valid_discard_ranked else (None, None)
    )
    best_keep, keep_note = (
        pick_best_with_preference(valid_keep_ranked)
        if valid_keep_ranked else (None, None)
    )

    # Attach per-entry notes so each row in the TSV only gets its own note
    best_discard_entry = {**best_discard, 'note': discard_note or ''} if best_discard else None
    best_keep_entry    = {**best_keep,    'note': keep_note    or ''} if best_keep    else None

    # ── apply selection logic ─────────────────────────────────────────────────
    selected = []
    warnings = []

    discard_above = best_discard and best_discard['mean_retained'] >= RETENTION_THRESHOLD
    keep_above    = best_keep    and best_keep['mean_retained']    >= RETENTION_THRESHOLD

    if discard_above and keep_above:
        # Both groups adequate: one from each
        selected.append({**best_discard_entry, 'group': 'best_discard'})
        selected.append({**best_keep_entry,    'group': 'best_keep'})

    elif discard_above:
        # Only discard adequate: up to two best discard above threshold
        above_threshold = [e for e in valid_discard_ranked
                           if e['mean_retained'] >= RETENTION_THRESHOLD]
        selected.extend(
            [{**e, 'group': f'best_discard_{i+1}',
            'note': discard_note if e['combo'] == best_discard_entry['combo'] else ''}
            for i, e in enumerate(above_threshold[:2])]
        )
        if not valid_keep:
            warnings.append(
                "No valid keep combos found (all had 100% retention with 0% "
                "adapter hits or 0% retention) — using discard combo(s) only."
            )
        elif not keep_above:
            warnings.append(
                f"Best keep combo '{best_keep['combo']}' had only "
                f"{best_keep['mean_retained']}% mean retention — "
                "using discard combo(s) only."
            )

    elif keep_above:
        # Only keep adequate: up to two best keep above threshold
        above_threshold = [e for e in valid_keep_ranked
                           if e['mean_retained'] >= RETENTION_THRESHOLD]
        selected.extend(
            [{**e, 'group': f'best_keep_{i+1}',
            'note': keep_note if e['combo'] == best_keep_entry['combo'] else ''}
             for i, e in enumerate(above_threshold[:2])]
        )
        if not valid_discard:
            warnings.append(
                "No valid discard combos found (all had 0% retention) — "
                "using keep combo(s) only."
            )
        elif not discard_above:
            warnings.append(
                f"Best discard combo '{best_discard['combo']}' had only "
                f"{best_discard['mean_retained']}% mean retention — "
                "falling back to keep combo(s)."
            )

    else:
        # Both below threshold: top two from combined valid pool
        warnings.append(
            "Both discard and keep combos are below the retention threshold "
            f"({RETENTION_THRESHOLD}%) — picking top two from combined valid pool. "
            "Primers may not be effectively trimmed."
        )
        combined = sorted(valid_discard + valid_keep, key=rank_combined, reverse=True)
        if not combined:
            # Absolute fallback — nothing valid at all, take raw best from each group
            combined = sorted(discard_combos + keep_combos,
                              key=rank_combined, reverse=True)
            warnings.append(
                "No valid combos found in either group — using best available."
            )
        selected.extend(
            [{**e, 'group': f'fallback_{i+1}', 'note': ''} for i, e in
             enumerate(combined[:2])]
        )

    # ── print summary ─────────────────────────────────────────────────────────
    for entry in selected:
        print(f"{entry['group']}: {entry['combo']} "
              f"(mean_retained={entry['mean_retained']}%, "
              f"mean_max_adapter={entry['mean_max_adapter']}%)")

    for w in warnings:
        print(f"Warning: {w}", file=sys.stderr)

    # ── write outputs ─────────────────────────────────────────────────────────
    out_cols = ['group', 'combo', 'mean_retained', 'mean_max_adapter', 'warning', 'note']

    with open(args.output_summary, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter='\t')
        w.writeheader()
        for entry in selected:
            w.writerow({
                'group':                  entry['group'],
                'combo':                  entry['combo'],
                'mean_retained':          entry['mean_retained'],
                'mean_max_adapter':       entry['mean_max_adapter'],
                'warning':                "; ".join(warnings) if warnings else "",
                'note':                   entry.get('note', ''),
            })

    # Resolve selected combo labels by group for Nextflow channel filtering
    discard_labels = [e['combo'] for e in selected if 'discard' in e['group']]
    keep_labels    = [e['combo'] for e in selected if 'keep'    in e['group']]

    # Fallback combos go to keep file if no explicit group assignment
    fallback_labels = [e['combo'] for e in selected if 'fallback' in e['group']]
    keep_labels    += fallback_labels

    with open(args.output_best_discard, 'w') as fh:
        fh.write('\n'.join(discard_labels) + '\n' if discard_labels else '')

    with open(args.output_best_keep, 'w') as fh:
        fh.write('\n'.join(keep_labels) + '\n' if keep_labels else '')

    print(f"Written: {args.output_summary}, "
          f"{args.output_best_discard}, {args.output_best_keep}")


if __name__ == '__main__':
    main()