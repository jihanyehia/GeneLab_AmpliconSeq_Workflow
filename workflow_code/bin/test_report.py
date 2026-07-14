#!/usr/bin/env python3
"""
Collect test-mode TSVs staged into the current working directory and render
a self-contained HTML report.

Usage:
    test_report.py --test-level <level> --f-primer <seq> --r-primer <seq>
                   --target-region <16S|18S|ITS> --input <accession or file>
                   --multiqc-json <path> --is-paired <true|false>
                   --auto-lt <int> --auto-rt <int>
"""

import argparse
import csv, os, glob, sys, json
from pathlib import Path
from datetime import datetime


def load_tsv(path):
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            rows = list(reader)
            return (reader.fieldnames or []), rows
    except Exception:
        return [], []


def collect_tsvs(pattern):
    cols, all_rows = [], []
    for f in sorted(glob.glob(pattern)):
        c, rows = load_tsv(f)
        if not cols and c:
            cols = c
        all_rows.extend(rows)
    return cols, all_rows


def rows_to_html(cols, rows,
                 warn_col = None, warn_threshold = 50.0,
                 selected_col = None, selected_vals = None):
    """
    Render a list of dicts as an HTML table.

    Row colouring priority (highest to lowest):
        selected  — row is one of the best combos selected for filter grid (green)
        warn      — value in warn_col is between 0 and warn_threshold (amber)
        fail      — value in warn_col is the string 'ERROR' (red)
    """
    if not rows:
        return "<p>No data collected.</p>"

    selected_vals = selected_vals or set()
    header_html = "".join(f"<th>{c}</th>" for c in cols)
    body_html = ""

    for row in rows:
        css_class = ""

        if selected_col and row.get(selected_col) in selected_vals:
            css_class = "selected"
        elif warn_col:
            val = row.get(warn_col, "")
            if val == "ERROR":
                css_class = "fail"
            elif val not in ("", None):
                try:
                    if float(val) > 0 and float(val) < warn_threshold:
                        css_class = "warn"
                except ValueError:
                    pass

        cells = "".join(f"<td>{row.get(c, '')}</td>" for c in cols)

        body_html += f'<tr class="{css_class}">{cells}</tr>'

    return f"""
    <table>
      <thead><tr>{header_html}</tr></thead>
      <tbody>{body_html}</tbody>
    </table>
    """


def build_quality_chart_html(multiqc_json_path, is_paired, auto_lt, auto_rt):
    """
    Parse MultiQC per-base sequence quality data and render interactive
    Plotly charts with per-sample quality curves, median line, Q20 threshold,
    and vertical lines at auto-picked truncLen positions.
    Returns an HTML string containing the charts, or a placeholder if data
    is unavailable.
    """
    try:
        with open(multiqc_json_path) as f:
            d = json.load(f)
        lines = (d["report_plot_data"]
                  ["fastqc_per_base_sequence_quality_plot"]
                  ["datasets"][0]
                  ["lines"])
    except Exception as e:
        return f"<p>Quality profile data unavailable: {e}</p>"

    def compute_median(dir_lines):
        """Compute median quality per position across all samples."""
        all_positions = sorted(set(p[0] for l in dir_lines for p in l['pairs']))
        median_ys = []
        for pos in all_positions:
            vals = sorted([p[1] for l in dir_lines for p in l['pairs'] if p[0] == pos])
            if vals:
                n = len(vals)
                median = (vals[n // 2] if n % 2 == 1
                          else (vals[n // 2 - 1] + vals[n // 2]) / 2)
                median_ys.append(median)
            else:
                median_ys.append(None)
        return all_positions, median_ys

    def build_traces(dir_lines, auto_trunc):
        traces = []
        all_positions, median_ys = compute_median(dir_lines)

        # Per-sample traces (light grey, hidden from legend)
        for line in dir_lines:
            xs = [p[0] for p in line['pairs']]
            ys = [p[1] for p in line['pairs']]
            traces.append({
                "x": xs, "y": ys,
                "type": "scatter", "mode": "lines",
                "name": line['name'],
                "line": {"color": "#d3d1c7", "width": 1},
                "opacity": 0.7,
                "showlegend": False,
                "hovertemplate": (
                    f"{line['name']}<br>"
                    "Position: %{x} bp<br>"
                    "Quality: %{y:.1f}<extra></extra>"
                )
            })

        # Median trace (purple)
        traces.append({
            "x": all_positions, "y": median_ys,
            "type": "scatter", "mode": "lines",
            "name": "Median",
            "line": {"color": "#534ab7", "width": 2},
            "hovertemplate": (
                "Median<br>"
                "Position: %{x} bp<br>"
                "Quality: %{y:.1f}<extra></extra>"
            )
        })

        # Q20 horizontal reference line (amber dashed)
        traces.append({
            "x": [min(all_positions), max(all_positions)],
            "y": [20, 20],
            "type": "scatter", "mode": "lines",
            "name": "Q20 threshold",
            "line": {"color": "#ef9f27", "width": 1, "dash": "dash"},
            "hoverinfo": "skip"
        })

        # Auto truncLen vertical line (coral dashed)
        if auto_trunc and auto_trunc > 0:
            traces.append({
                "x": [auto_trunc, auto_trunc],
                "y": [0, 40],
                "type": "scatter", "mode": "lines",
                "name": f"Auto truncLen ({auto_trunc})",
                "line": {"color": "#d85a30", "width": 1.5, "dash": "dash"},
                "hoverinfo": "skip"
            })

        return traces, all_positions

    def build_chart(dir_lines, direction, auto_trunc, div_id):
        if not dir_lines:
            return f"<p>No {direction} quality data found.</p>"

        traces, all_positions = build_traces(dir_lines, auto_trunc)

        layout = {
            "title": {
                "text": direction,
                "font": {"size": 13, "color": "#3c3489"}
            },
            "xaxis": {
                "title": "Position (bp)",
                "range": [0, max(all_positions)],
                "gridcolor": "#f1efe8",
                "linecolor": "#d3d1c7",
                "tickfont": {"size": 11}
            },
            "yaxis": {
                "title": "Phred score",
                "range": [0, 40],
                "gridcolor": "#f1efe8",
                "linecolor": "#d3d1c7",
                "tickfont": {"size": 11}
            },
            "plot_bgcolor":  "#fff",
            "paper_bgcolor": "#fff",
            "legend": {"font": {"size": 11}},
            "margin": {"l": 50, "r": 20, "t": 36, "b": 50},
            "hovermode": "x unified"
        }

        traces_json = json.dumps(traces)
        layout_json = json.dumps(layout)

        return f"""
        <div id="{div_id}"
             style="width:100%;height:300px;border:1px solid #d3d1c7;
                    border-radius:8px;overflow:hidden;"></div>
        <script>
        Plotly.newPlot("{div_id}", {traces_json}, {layout_json},
            {{responsive: true, displayModeBar: false}});
        </script>
        """

    r1_lines = [l for l in lines if l['name'].endswith('_R1')]
    r2_lines = [l for l in lines if l['name'].endswith('_R2')] if is_paired else []

    html  = build_chart(r1_lines, "R1", auto_lt, "quality_chart_r1")
    if is_paired and r2_lines:
        html += build_chart(r2_lines, "R2", auto_rt, "quality_chart_r2")

    return html


def build_filter_html(filter_cols, filter_rows):
    """
    Build the filter read-count tracking HTML, grouped by cutadapt_label so
    discard vs keep results are shown in separate tables.
    If cutadapt was skipped, no label is found and no header is rendered.
    truncLen_source is rendered as a badge in the first column.
    Rows are visually grouped by parameter combo with alternating backgrounds.
    """
    if not filter_rows:
        return "<p>Not run at this test level.</p>"

    combo_cols = ['left_trunc', 'right_trunc', 'left_maxEE', 'right_maxEE']
    
    # Exclude cutadapt_label and truncLen_source from regular columns:
    # cutadapt_label is used for grouping, truncLen_source is rendered as a badge
    display_cols = [c for c in filter_cols if c not in ('cutadapt_label', 'truncLen_source')]

    source_labels = {
        'default': ('Default', 'badge-default'),
        'auto':    ('Auto',    'badge-auto'),
        'user':    ('User',    'badge-user'),
        'both':    ('Both',    'badge-both'),
    }

    # Group rows by cutadapt_label
    groups = {}
    for row in filter_rows:
        label = row.get('cutadapt_label', '')
        groups.setdefault(label, []).append(row)

    # Determine whether cutadapt was actually run:
    # if there is only one group and its label is empty or a placeholder, skip sub-headers
    unique_labels = list(groups.keys())
    show_subheaders = not (
        len(unique_labels) == 1 and unique_labels[0] in ('', 'none', 'raw', None)
    )

    html = ""
    for label, rows in sorted(groups.items()):
        if show_subheaders:
            html += f"<h3>{label}</h3>"
        
        # truncLen_source badge is the first column
        header_html = (
            "<th>truncLen source</th>"
            + "".join(f"<th>{c}</th>" for c in display_cols)
        )
        body_html = ""
        combo_index = 0
        last_combo = None

        for row in rows:
            current_combo = tuple(row.get(c, '') for c in combo_cols)
            if current_combo != last_combo:
                if last_combo is not None:
                    combo_index += 1
                last_combo = current_combo

            # Alternate row background for each unique combo
            css_class = "combo-alt" if combo_index % 2 == 1 else ""

            # Override with fail highlight if this row has an error note —
            # takes priority over combo alternation since failures need immediate attention
            if row.get('notes', '').strip():
                css_class = "fail"

            src = row.get('truncLen_source', '')
            src_label, src_css = source_labels.get(src, (src, 'badge-neutral'))
            badge_cell = f'<td><span class="badge {src_css}">{src_label}</span></td>'

            cells = badge_cell + "".join(f"<td>{row.get(c, '')}</td>" for c in display_cols)
            body_html += f'<tr class="{css_class}">{cells}</tr>'

        html += f"""
        <table>
          <thead><tr>{header_html}</tr></thead>
          <tbody>{body_html}</tbody>
        </table>
        """

    return html


def build_recommendations(cutadapt_rows, best_combos):
    """
    Build an HTML recommendations list from the best cutadapt combos.
    """
    recs = []

    # Use pre-selected best combos from TEST_PICK_BEST_CUTADAPT rather than re-deriving them here
    if best_combos:
        for row in best_combos:
            group = row.get('group', '')
            combo = row.get('combo', '')
            ret   = row.get('mean_retained', '')
            adp   = row.get('mean_max_adapter', '')
            label = 'discard untrimmed' if group == 'best_discard' else 'keep untrimmed'
            recs.append(
                f"Best cutadapt combo ({label}): "
                f"<code>{combo}</code> — "
                f"{ret}% mean reads retained, "
                f"{adp}% mean reads with adapter."
            )

    recs_html = ("".join(f"<li>{r}</li>" for r in recs)
                if recs else "<li>No recommendations available at this test level.</li>")

    return recs_html


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--test-level", required=True)
    parser.add_argument("--f-primer", required=True)
    parser.add_argument("--r-primer", required=True)
    parser.add_argument("--target-region", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--multiqc-json",  required=True,
                        help="Path to multiqc_data.json from raw MultiQC")
    parser.add_argument("--is-paired",     required=True,
                        help="Whether data is paired-end (true/false)")
    parser.add_argument("--auto-lt",       required=False, type=int, default=None,
                        help="Auto-picked left truncLen from PICK_TRUNCLEN")
    parser.add_argument("--auto-rt",       required=False, type=int, default=None,
                        help="Auto-picked right truncLen from PICK_TRUNCLEN")

    args = parser.parse_args()

    # ── load data ─────────────────────────────────────────────────────────────────
    qcols, qrows             = load_tsv("sample_quality_rank.tsv")
    pcols, prows             = load_tsv("primer_detection_summary.tsv")
    cut_cols, cut_rows       = collect_tsvs("*_cutadapt_grid.tsv")
    best_cols, best_combos   = load_tsv("best_cutadapt_combos.tsv")
    filter_cols, filter_rows = collect_tsvs("filter_grid_*_read_counts_tracking.tsv")

    # Resolve best combo labels for highlighting in the cutadapt table
    best_combo_labels = {row.get('combo', '') for row in best_combos}

    primer_decision = "not run"
    try:
        with open("primer_decision.txt") as fh:
            primer_decision = fh.read().strip()
    except FileNotFoundError:
        pass

    # ── render sections ───────────────────────────────────────────────────────────
    if qrows:
        quality_html = rows_to_html(qcols, qrows,
                                   selected_col='selected', selected_vals={'yes'})
        quality_note = (
            "Samples were ranked by composite MultiQC quality metrics and selected "
            "from the provided input dataset."
        )
    else:
        quality_html = "<p>No quality ranking data available.</p>"
        quality_note = (
            "Samples were provided directly via the input file for test mode; "
            "no additional sample selection was performed."
        )

    primer_html = (rows_to_html(pcols, prows)
                if prows else "<p>Not run at this test level.</p>")

    if cut_rows:
        # Highlight best combos in green, warn on low retention in amber
        cutadapt_html = rows_to_html(cut_cols, cut_rows, 
                                     warn_col='pct_retained', warn_threshold=50.0,
                                     selected_col='combo', selected_vals=best_combo_labels)
        if best_combos:
            # Extract selected combo labels, both can have multiple entries
            discard_combos_selected = [r['combo'] for r in best_combos if 'discard' in r.get('group', '')]
            keep_combos_selected    = [r['combo'] for r in best_combos if 'keep'    in r.get('group', '')]

            discard_label = ", ".join(f"<code>{c}</code>" for c in discard_combos_selected) \
                    if discard_combos_selected else "<em>not selected</em>"
            keep_label    = ", ".join(f"<code>{c}</code>" for c in keep_combos_selected) \
                            if keep_combos_selected else "<em>not selected</em>"
            
            # Collect all warning messages, if any, from the best_combos rows
            warning_msgs = [r['warning'] for r in best_combos
                    if r.get('warning', '').strip()]
            warning_html = "".join(
                f'<p class="warn-note">&#9888; {msg}</p>'
                for msg in warning_msgs
            )

            # Collect preference tiebreaker notes from best_combos rows
            note_msgs = [r['note'] for r in best_combos if r.get('note', '').strip()]
            note_html = "".join(
                f'<p class="pref-note">&#8505; {msg}</p>'
                for msg in note_msgs
            )

            cutadapt_html += f"""
            <p class="note">
                <span class="badge badge-ok">&#10003; Selected for filter grid</span>
                &nbsp;
                <b>Discard:</b> {discard_label}
                &nbsp;&nbsp;
                <b>Keep:</b> {keep_label}
            </p>
            {warning_html}
            {note_html}
            """
    else:
        if args.test_level in ('primers', 'full') and primer_decision == 'primers_not_found':
            cutadapt_html = (
                "<p>No primers were detected in raw reads; cutadapt was skipped.</p>"
            )
        else:
            cutadapt_html = "<p>Not run at this test level.</p>"

    # Quality profile charts only when the filter grid runs (truncLen is only relevant in that context)
    if args.test_level in ('filter', 'full'):
        quality_chart_html = build_quality_chart_html(
            args.multiqc_json, args.is_paired, args.auto_lt, args.auto_rt
        )
    else:
        quality_chart_html = "<p>Not run at this test level.</p>"

    filter_html = build_filter_html(filter_cols, filter_rows)

    its_note = ""
    if args.target_region == "ITS":
        its_note = """
        <p class="note" style="background:#faeeda;border:1px solid #ef9f27;border-radius:6px;padding:8px 12px;margin-top:8px;">
            &#9888; ITS amplicons have highly variable lengths due to the nature of the target region.
            Truncating reads to a fixed length is generally not recommended for ITS data as it can
            disproportionately discard reads from taxa with longer ITS sequences.
            The <span class="badge badge-default">Default</span> (0/0, no truncation) combo is
            typically the most appropriate choice for ITS.
        </p>
        """

    decision_badge_cls = {
        'primers_found':     'badge-ok',
        'primers_not_found': 'badge-warn',
        'not run':           'badge-neutral',
    }.get(primer_decision, 'badge-neutral')
    
    decision_label = primer_decision.replace('_', ' ')

    recommendations_html = build_recommendations(cut_rows, best_combos)

    # ── render HTML ───────────────────────────────────────────────────────────────
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    html = f"""\
    <!DOCTYPE html>
    <html lang="en">
    <head>
    <meta charset="utf-8">
    <title>AmpIllumina Test Mode Report</title>
    <script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.34.0/plotly.min.js"></script>
    <style>
    body      {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
                margin: 0; padding: 20px 40px; background: #f8f8f6;
                color: #2c2c2a; font-size: 14px; }}
    h1        {{ font-size: 22px; font-weight: 500;
                border-bottom: 1px solid #d3d1c7; padding-bottom: 8px; }}
    h2        {{ font-size: 17px; font-weight: 500; margin-top: 32px; color: #3c3489; }}
    h3        {{ font-size: 14px; font-weight: 500; margin-top: 20px; color: #5f5e5a; }}
    .meta     {{ background: #fff; border: 1px solid #d3d1c7; border-radius: 8px;
                padding: 12px 20px; margin-bottom: 24px;
                display: flex; gap: 32px; flex-wrap: wrap; }}
    .meta span {{ font-size: 13px; color: #5f5e5a; }}
    .meta b    {{ color: #2c2c2a; }}
    table     {{ border-collapse: collapse; width: 100%; background: #fff;
                border-radius: 8px; overflow: hidden; border: 1px solid #d3d1c7;
                margin-top: 8px; font-size: 13px; }}
    th        {{ background: #eeedfe; color: #3c3489; font-weight: 500;
                text-align: left; padding: 8px 12px;
                border-bottom: 1px solid #afa9ec; }}
    td        {{ padding: 6px 12px; border-bottom: 1px solid #f1efe8; }}
    tr:last-child td  {{ border-bottom: none; }}
    tr.fail td        {{ background: #fcebeb; color: #a32d2d; }}
    tr.warn td        {{ background: #faeeda; color: #854f0b; }}
    tr.selected td    {{ background: #eaf3de; color: #3b6d11; font-weight: 500; }}
    tr.combo-alt td   {{ background: #f5f4f0; }}
    .badge    {{ display: inline-block; padding: 2px 10px; border-radius: 12px;
                font-size: 12px; font-weight: 500; border: 1px solid transparent; }}
    .badge-ok      {{ background: #eaf3de; color: #3b6d11; border-color: #97c459; }}
    .badge-warn    {{ background: #faeeda; color: #854f0b; border-color: #ef9f27; }}
    .badge-neutral {{ background: #f1efe8; color: #5f5e5a; border-color: #b4b2a9; }}
    .badge-default {{ background: #f1efe8; color: #444441; border-color: #b4b2a9; }}
    .badge-auto    {{ background: #e6f1fb; color: #0c447c; border-color: #85b7eb; }}
    .badge-user    {{ background: #eeedfe; color: #3c3489; border-color: #afa9ec; }}
    .badge-both    {{ background: #eaf3de; color: #3b6d11; border-color: #97c459; }}
    dl.badge-legend {{ display: grid; grid-template-columns: auto 1fr;
                        gap: 4px 12px; margin: 8px 0 0; }}
    dl.badge-legend dt {{ display: flex; align-items: center; }}
    dl.badge-legend dd {{ font-size: 13px; color: #5f5e5a; margin: 0;
                        display: flex; align-items: center; }}
    .recs     {{ background: #e6f1fb; border: 1px solid #85b7eb; border-radius: 8px;
                padding: 12px 20px; margin-top: 8px; }}
    .recs li  {{ margin: 4px 0; color: #0c447c; font-size: 13px; }}
    code      {{ background: #f1efe8; padding: 1px 6px; border-radius: 4px;
                font-family: 'SF Mono', 'Consolas', monospace; font-size: 12px; }}
    section   {{ margin-bottom: 28px; }}
    p.note    {{ color: #5f5e5a; font-size: 13px; }}
    p.warn-note    {{ color: #854f0b; font-size: 13px; }}
    p.pref-note {{ background: #e6f1fb; border: 1px solid #85b7eb; border-radius: 6px;
              padding: 8px 12px; color: #0c447c; font-size: 13px; margin-top: 6px; }}
    .chart-grid    {{ display: grid;
                    grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
                    gap: 16px; margin-top: 8px; }}
    </style>
    </head>
    <body>
    <h1>AmpIllumina — Test Mode Report</h1>

    <div class="meta">
    <span><b>Generated:</b> {now}</span>
    <span><b>Target region:</b> {args.target_region}</span>
    <span><b>Input:</b> {args.input}</span>
    <span><b>Test level:</b> {args.test_level}</span>
    <span><b>F primer:</b> <code>{args.f_primer}</code></span>
    <span><b>R primer:</b> <code>{args.r_primer}</code></span>
    <span><b>Primer detection:</b>
        <span class="badge {decision_badge_cls}">{decision_label}</span>
    </span>
    </div>

    <section>
    <h2>1. Sample selection</h2>
    <p class="note">
        {quality_note}
    </p>
    {quality_html}
    </section>

    <section>
    <h2>2. Primer detection in raw reads</h2>
    <p class="note">
        Searches the first N reads of each selected sample for the forward and
        reverse primers and their reverse complements. 
        High hit rates in any orientation confirm primers are present and indicate 
        the expected trimming strategy. 
        High F_primer_pct on R1 and R_primer_rc_pct on R2 confirms correct primer 
        orientation.
    </p>
    {primer_html}
    </section>

    <section>
    <h2>3. Cutadapt parameter grid</h2>
    <p class="note">
        All 8 combinations of anchored/unanchored &times; linked/unlinked &times;
        discard_untrimmed, tested on each selected sample.
        Rows highlighted <span style="color:#3b6d11;font-weight:500">green</span>
        were selected as the best discard and best keep combo for the filter grid.
        Rows with &lt;50% read retention are highlighted amber.
    </p>
    {cutadapt_html}
    </section>

    <section>
    <h2>4. Raw read quality profiles</h2>
    <p class="note">
        Per-base mean quality scores across all samples from raw MultiQC data.
        Each grey line represents one sample; the purple line is the median across
        all samples. The amber dashed line marks the Q20 threshold used to
        auto-pick truncLen. The coral dashed line marks the auto-picked truncLen
        position applied in the filter grid below.
        Hover over the chart to inspect exact quality values per sample and position.
    </p>
    <p class="note" style="background:#faeeda;border:1px solid #ef9f27;border-radius:6px;padding:8px 12px;margin-top:4px;">
        &#9888; truncLen values are auto-picked from <b>raw</b> read quality profiles.
        After primer trimming, reads may be shorter than the auto-picked truncLen,
        particularly for R2 reads where primers are trimmed from the 5' end.
        If auto truncLen combos show poor retention in the filter grid below, consider
        using the <span class="badge badge-default">Default</span> (0/0) combo or
        supplying manual values via <code>--test_trunc_left</code> / <code>--test_trunc_right</code>.
    </p>
    <div class="chart-grid">
    {quality_chart_html}
    </div>
    </section>

    <section>
    <h2>5. DADA2 filtering parameter grid</h2>
    <p class="note">
        Per-sample DADA2 counts tracking across all truncLen &times; maxEE
        combinations run on reads trimmed by the best discard and best keep
        cutadapt combos (or raw reads when applicable).
        The <b>truncLen source</b> badge indicates how each truncLen value was chosen:
    </p>
    <dl class="badge-legend">
        <dt><span class="badge badge-default">Default</span></dt>
        <dd>no truncation baseline (0/0, always included)</dd>
        <dt><span class="badge badge-auto">Auto</span></dt>
        <dd>picked from raw MultiQC Q20 quality threshold</dd>
        <dt><span class="badge badge-user">User</span></dt>
        <dd>supplied via <code>--test_trunc_left</code> / <code>--test_trunc_right</code></dd>
        <dt><span class="badge badge-both">Both</span></dt>
        <dd>matches both auto and user</dd>
    </dl>
    {its_note}
    {filter_html}
    </section>

    <section>
    <h2>Recommendations</h2>
    <div class="recs"><ul>{recommendations_html}</ul></div>
    </section>
    </body>
    </html>"""

    outfile = f"{args.input}_{args.target_region}_test_mode_report.html"
    with open(outfile, "w") as fh:
        fh.write(html)

    print(f"Written: {outfile}")


if __name__ == "__main__":
    main()