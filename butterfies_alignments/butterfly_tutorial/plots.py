"""Plot measured BioseqDB query results without running an external aligner.

All interval coordinates are zero-based and half-open, relative to the imported
sequence. Database access belongs in workflow.py; these functions only draw the
records returned by that workflow. No observations are synthesized here.
"""
from __future__ import annotations

import html
import json
import math
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

import plotly.graph_objects as go
from plotly.subplots import make_subplots


PALETTES = {
    "dark": {
        "background": "#071320", "panel": "#0a1826", "text": "#eef4f9",
        "muted": "#a9bccb", "grid": "#20364a", "cyan": "#65b9dd",
        "mint": "#74c8a8", "orange": "#f4894f", "purple": "#c5a8d6",
        "red": "#ef7a5a", "yellow": "#e2c478",
    },
    "light": {
        "background": "#f5f8fc", "panel": "#ffffff", "text": "#173247",
        "muted": "#405c72", "grid": "#ccd9e4", "cyan": "#146c91",
        "mint": "#167251", "orange": "#b44d1c", "purple": "#7d468c",
        "red": "#b3422b", "yellow": "#805b16",
    },
}
BASE_COLORS = {"A": "purple", "C": "cyan", "G": "orange", "T": "mint", "N": "muted"}
Record = Mapping[str, Any]


def records(values: Iterable[Record]) -> list[dict[str, Any]]:
    """Accept ordinary records or a pandas DataFrame without requiring pandas."""
    if hasattr(values, "to_dict"):
        values = values.to_dict(orient="records")
    return [dict(row) for row in values]


def _number(value: Any, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _boolean(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value in (0, 1):
        return bool(value)
    raise ValueError("is_reverse must be a Boolean, not a text value")


def validate_inputs(sequences: Iterable[Record], matches: Iterable[Record],
                    gc_windows: Iterable[Record], genes: Iterable[Record] = ()):
    """Reject empty or incompatible inputs before a misleading chart is written."""
    seq, hits, windows, annotations = map(records, (sequences, matches, gc_windows, genes))
    if len(seq) < 2:
        raise ValueError("At least two imported sequences are required")
    by_id = {str(row["accession"]): row for row in seq}
    if len(by_id) != len(seq):
        raise ValueError("Sequence accession values must be unique")
    for row in seq:
        length = int(row["length"])
        if length <= 0 or length != row["length"]:
            raise ValueError("Sequence lengths must be positive integers")
        counts = [int(row[base.lower()]) for base in "ACGTN"]
        if any(n < 0 for n in counts) or sum(counts) != length:
            raise ValueError(f"ACGTN counts do not sum to {row['accession']} length")
        called = sum(counts[:4])
        gc = row.get("gc_percent")
        if called and (gc is None or abs(_number(gc, "GC") - 100 * (counts[1] + counts[2]) / called) > .011):
            raise ValueError("GC must use A+C+G+T as its denominator")
        if not called and gc is not None:
            raise ValueError("An all-N sequence must have null GC")

    def check_interval(row: Record, accession: str, start: str, end: str):
        if accession not in by_id:
            raise ValueError(f"Unknown accession {accession}")
        lo, hi = int(row[start]), int(row[end])
        if lo != row[start] or hi != row[end] or not 0 <= lo < hi <= int(by_id[accession]["length"]):
            raise ValueError(f"Out-of-bounds interval for {accession}: [{lo}, {hi})")

    if not hits:
        raise ValueError("BioseqDB returned no alignment rows; charts were not generated")
    pairs = set()
    for row in hits:
        target, query = str(row["target_accession"]), str(row["query_accession"])
        pairs.add((target, query))
        check_interval(row, target, "target_start", "target_end")
        check_interval(row, query, "query_start", "query_end")
        row["is_reverse"] = _boolean(row["is_reverse"])
        if row.get("identity_percent") is not None and not 0 <= _number(row["identity_percent"], "identity") <= 100:
            raise ValueError("Identity must be between zero and 100 percent")
    if len(pairs) != 1:
        raise ValueError("The pairwise figure expects one target/query accession pair")
    if not windows:
        raise ValueError("No GC windows were returned")
    for row in windows:
        check_interval(row, str(row["accession"]), "start", "end")
        if row.get("gc_percent") is not None and not 0 <= _number(row["gc_percent"], "window GC") <= 100:
            raise ValueError("Window GC must be between zero and 100 percent")
        if not 0 <= _number(row["n_fraction"], "unknown-base fraction") <= 1:
            raise ValueError("n_fraction must be in [0, 1]")
    for row in annotations:
        check_interval(row, str(row["accession"]), "start", "end")
    return seq, hits, windows, annotations


def _layout(fig: go.Figure, title: str, theme: str = "dark", *, height: int = 520,
            subtitle: str = "", left: int = 84) -> go.Figure:
    p = PALETTES[theme]
    fig.update_layout(
        template="none", width=860, height=height,
        title={"text": html.escape(title) + (f"<br><sup>{html.escape(subtitle)}</sup>" if subtitle else ""),
               "font": {"size": 23}, "x": .025, "xanchor": "left"},
        font={"family": "Arial, sans-serif", "size": 17, "color": p["text"]},
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        margin={"t": 116 if subtitle else 88, "b": 80, "l": left, "r": 32},
        legend={"orientation": "h", "y": -0.2, "x": 0, "font": {"size": 15}},
        hoverlabel={"bgcolor": p["panel"], "font": {"color": p["text"], "size": 14}},
        modebar={"color": p["muted"], "activecolor": p["text"], "bgcolor": p["panel"]},
        hovermode="closest",
    )
    fig.update_xaxes(gridcolor=p["grid"], zeroline=False, linecolor=p["grid"],
                     tickfont={"color": p["muted"], "size": 14}, title_font={"size": 16}, automargin=True)
    fig.update_yaxes(gridcolor=p["grid"], zeroline=False, linecolor=p["grid"],
                     tickfont={"color": p["muted"], "size": 14}, title_font={"size": 16}, automargin=True)
    return fig


def _label(row: Record) -> str:
    organism = str(row.get("organism", row["accession"]))
    if organism.startswith("Heliconius "):
        organism = "H. " + organism.removeprefix("Heliconius ")
    return f"{html.escape(organism)}<br>{html.escape(str(row['accession']))}"


def _sequence_color(row: Record) -> str:
    return "mint" if "erato" in str(row.get("organism", "")).casefold() else "cyan"


def sequence_lengths(seq: list[Record], theme="dark") -> go.Figure:
    p = PALETTES[theme]
    fig = go.Figure(go.Bar(
        x=[r["length"] for r in seq], y=[_label(r) for r in seq], orientation="h",
        marker_color=[p[_sequence_color(r)] for r in seq],
        text=[f"{int(r['length']):,} bp" for r in seq], textposition="outside", cliponaxis=False,
        hovertemplate="%{y}<br>Imported length: %{x:,} bp<extra></extra>"))
    _layout(fig, "The DNA intervals we imported", theme, height=max(410, 180 + 95 * len(seq)),
            subtitle="Lengths describe selected intervals, not complete species genomes", left=240)
    fig.update_xaxes(title="Sequence length (bp)", range=[0, max(r["length"] for r in seq) * 1.29])
    fig.update_yaxes(autorange="reversed", showgrid=False)
    return fig


def base_counts(seq: list[Record], theme="dark") -> go.Figure:
    p = PALETTES[theme]
    fig = go.Figure()
    for base, color in BASE_COLORS.items():
        fig.add_bar(name=base, x=[_label(r) for r in seq], y=[r[base.lower()] for r in seq],
                    marker_color=p[color], customdata=[[100 * r[base.lower()] / r["length"]] for r in seq],
                    hovertemplate=f"{base}: %{{y:,}} bases<br>%{{customdata[0]:.2f}}% of interval<extra></extra>")
    _layout(fig, "How many of each DNA letter?", theme, height=540,
            subtitle="A, C, G and T are called bases; N marks an unresolved position")
    fig.update_layout(barmode="stack", legend={"y": -0.27})
    fig.update_yaxes(title="Number of bases")
    fig.update_xaxes(showgrid=False, tickfont={"size": 13})
    return fig


def source_record_lengths(seq: list[Record], theme="dark") -> go.Figure:
    p = PALETTES[theme]
    fig = go.Figure()
    for field, label, color in (("source_length", "Complete NCBI record", "muted"), ("length", "Imported interval", "cyan")):
        fig.add_bar(x=[_label(r) for r in seq], y=[r[field] for r in seq], name=label,
            marker_color=p[color], customdata=[[r.get("source_role", "sequence record")] for r in seq],
            hovertemplate="%{y:,} bp<br>Record role: %{customdata[0]}<extra>%{fullData.name}</extra>")
    _layout(fig, "The selected windows sit inside larger records", theme, height=540,
            subtitle="One scaffold and one chromosome record; neither is a complete species genome")
    fig.update_layout(barmode="group", legend={"y": -.27})
    fig.update_yaxes(title="Record / interval length (bp)")
    fig.update_xaxes(showgrid=False, tickfont={"size": 13})
    return fig


def gc_content(seq: list[Record], theme="dark") -> go.Figure:
    p = PALETTES[theme]
    fig = go.Figure(go.Bar(x=[_label(r) for r in seq], y=[r["gc_percent"] for r in seq],
        marker_color=p["purple"], text=[f"{float(r['gc_percent']):.2f}%" if r["gc_percent"] is not None else "undefined" for r in seq],
        textposition="outside", customdata=[[r["n"], r["length"]] for r in seq],
        hovertemplate="GC: %{y:.3f}% of A/C/G/T<br>N bases: %{customdata[0]:,}<br>Total: %{customdata[1]:,}<extra></extra>"))
    _layout(fig, "GC content is a composition summary", theme, height=480,
            subtitle="100 × (G + C) / (A + C + G + T); unknown bases are excluded")
    fig.update_yaxes(title="GC among called bases (%)", range=[0, 100], ticksuffix="%")
    fig.update_xaxes(showgrid=False, tickfont={"size": 13})
    return fig


def composition_windows(seq: list[Record], windows: list[Record], theme="dark") -> go.Figure:
    p = PALETTES[theme]
    fig = make_subplots(rows=len(seq), cols=1, shared_xaxes=False,
        vertical_spacing=.17, subplot_titles=[_label(r) for r in seq],
        specs=[[{"secondary_y": True}] for _ in seq])
    for i, sequence in enumerate(seq, 1):
        rows = sorted((r for r in windows if r["accession"] == sequence["accession"]), key=lambda r: r["start"])
        if not rows:
            raise ValueError(f"No composition windows for {sequence['accession']}")
        x = [(r["start"] + r["end"]) / 2 for r in rows]
        custom = [[r["start"], r["end"]] for r in rows]
        fig.add_trace(go.Scatter(x=x, y=[r["gc_percent"] for r in rows], mode="lines", name="GC / called bases",
            connectgaps=False, line={"color": p[_sequence_color(sequence)], "width": 2}, customdata=custom, showlegend=i == 1,
            hovertemplate="[%{customdata[0]:,}, %{customdata[1]:,})<br>GC: %{y:.2f}%<extra></extra>"), row=i, col=1)
        fig.add_trace(go.Scatter(x=x, y=[100 * r["n_fraction"] for r in rows], mode="lines", name="Unknown bases (N)",
            line={"color": p["orange"], "width": 1.5, "dash": "dot"}, customdata=custom, showlegend=i == 1,
            hovertemplate="[%{customdata[0]:,}, %{customdata[1]:,})<br>N: %{y:.2f}%<extra></extra>"), row=i, col=1, secondary_y=True)
        fig.update_yaxes(range=[0, 100], title="GC (%)", row=i, col=1, secondary_y=False)
        fig.update_yaxes(range=[0, 100], title="N (%)", showgrid=False, row=i, col=1, secondary_y=True)
        fig.update_xaxes(title="Position in imported interval (bp)" if i == len(seq) else None, row=i, col=1)
    _layout(fig, "Composition changes along each interval", theme, height=260 + 270 * len(seq),
            subtitle="Each point summarizes one SQL window; all-N windows have no GC value")
    fig.update_annotations(font={"size": 15, "color": p["text"]})
    fig.update_layout(margin={"r": 90, "t": 158}, legend={"y": -0.12}, height=300 + 270 * len(seq))
    return fig


def _hit_text(hit: Record) -> str:
    out = [f"Target: [{int(hit['target_start']):,}, {int(hit['target_end']):,})",
           f"Query: [{int(hit['query_start']):,}, {int(hit['query_end']):,})",
           "Orientation: reverse" if hit["is_reverse"] else "Orientation: forward"]
    if hit.get("identity_percent") is not None:
        out.append(f"Verified identity: {float(hit['identity_percent']):.2f}%")
    if hit.get("score") is not None:
        out.append(f"Alignment score: {hit['score']}")
    return "<br>".join(html.escape(part) for part in out)


def dot_plot(hits: list[Record], theme="dark") -> go.Figure:
    p = PALETTES[theme]
    fig = go.Figure()
    for reverse, color, name in ((False, "cyan", "Forward"), (True, "orange", "Reverse")):
        selected = [r for r in hits if r["is_reverse"] == reverse]
        if not selected:
            continue
        x, y, hover = [], [], []
        for hit in selected:
            q0, q1 = hit["query_start"], hit["query_end"]
            if reverse:
                q0, q1 = q1, q0
            x += [hit["target_start"], hit["target_end"], None]
            y += [q0, q1, None]
            hover += [_hit_text(hit), _hit_text(hit), None]
        fig.add_scatter(x=x, y=y, mode="lines", name=name, text=hover,
                        line={"color": p[color], "width": 2, "dash": "dot" if reverse else "solid"},
                        hovertemplate="%{text}<extra></extra>")
    _layout(fig, "Where the intervals align", theme, height=650,
            subtitle=f"{len(hits):,} BioseqDB alignments; line direction indicates orientation")
    fig.update_xaxes(title=f"Target {hits[0]['target_accession']} position (bp)")
    fig.update_yaxes(title=f"Query {hits[0]['query_accession']} position (bp)")
    return fig


def _is_optix(gene: Record) -> bool:
    # NCBI may use a locus tag as the gene name; annotations.py verifies this
    # flag from the annotation metadata instead of guessing from plot position.
    return gene.get("is_optix") is True or str(gene["name"]).casefold() == "optix"


def _anchors(genes: list[Record], target: str, query: str) -> dict[str, float]:
    """Use one verified optix annotation per accession; refuse ambiguous anchors."""
    anchors = {}
    for accession in (target, query):
        candidates = [g for g in genes if g["accession"] == accession and _is_optix(g)]
        if len(candidates) != 1:
            raise ValueError(f"Need exactly one optix gene interval for {accession}, got {len(candidates)}")
        anchors[accession] = (candidates[0]["start"] + candidates[0]["end"]) / 2
    return anchors


def _gene_strand(value: Any) -> int:
    if value in ("+", 1, "+1", "1"):
        return 1
    if value in ("-", -1, "-1"):
        return -1
    raise ValueError("Gene strand must be + or - (or integer +1 / -1)")


def _display_transforms(genes: list[Record], target: str, query: str, centered: bool):
    """Return centers and signs; query reflection aligns optix gene orientation.

    x_target = target_position - target_center
    x_query = orientation * (query_position - query_center)

    orientation is -1 only when the two optix gene strands differ. This is
    equivalent to first mapping q to L-q and then subtracting L-query_center.
    It changes the display, not genomic coordinates or reported strand calls.
    """
    if not centered:
        return {target: 0, query: 0}, {target: 1, query: 1}
    centers = _anchors(genes, target, query)
    strands = {accession: _gene_strand(next(g["strand"] for g in genes
        if g["accession"] == accession and _is_optix(g)))
        for accession in (target, query)}
    return centers, {target: 1, query: 1 if strands[target] == strands[query] else -1}


def ribbons(seq: list[Record], hits: list[Record], genes: list[Record] = (), theme="dark", *,
            first_only: bool = False, centered: bool = False,
            window: tuple[float, float] | None = None) -> go.Figure:
    """Draw complete match polygons; axis clipping retains crossing alignments.

    Centered mode subtracts each species' optix midpoint independently and
    reflects the query axis when optix strands differ. It is a display
    transform, not a coordinate lift-over or an extra alignment.
    """
    p = PALETTES[theme]
    by_id = {row["accession"]: row for row in seq}
    target, query = hits[0]["target_accession"], hits[0]["query_accession"]
    shifts, signs = _display_transforms(genes, target, query, centered)
    def display(accession, position):
        return signs[accession] * (position - shifts[accession])
    def overlaps(hit, side, accession, lo, hi):
        a, b = sorted((display(accession, hit[f"{side}_start"]), display(accession, hit[f"{side}_end"])))
        return b > lo and a < hi
    selected = sorted(hits, key=lambda r: (-float(r.get("score") or 0), r["target_start"], r["query_start"]))
    if first_only:
        selected = selected[:1]
    if window:
        lo, hi = window
        selected = [h for h in selected if overlaps(h, "target", target, lo, hi) and overlaps(h, "query", query, lo, hi)]
        if not selected:
            raise ValueError("No alignment intersects both sides of the requested zoom window")
    fig = go.Figure()
    for reverse, color, name in ((False, "cyan", "Forward alignment"), (True, "orange", "Reverse alignment")):
        x, y, text = [], [], []
        for hit in selected:
            if hit["is_reverse"] != reverse:
                continue
            t0, t1 = display(target, hit["target_start"]), display(target, hit["target_end"])
            q0, q1 = display(query, hit["query_start"]), display(query, hit["query_end"])
            if reverse:
                q0, q1 = q1, q0
            x += [t0, t1, q1, q0, t0, None]
            y += [1, 1, 0, 0, 1, None]
            text += [_hit_text(hit)] * 5 + [None]
        if not x:
            continue
        fig.add_scatter(x=x, y=y, mode="lines", fill="toself", name=name,
            line={"color": p[color], "width": .9}, fillcolor=p[color],
            opacity=.7 if first_only else (.4 if theme == "dark" else .27), showlegend=True,
            text=text, hovertemplate="%{text}<extra></extra>")
    for accession, lane, color in ((target, 1, "cyan"), (query, 0, "mint")):
        extent = sorted((display(accession, 0), display(accession, by_id[accession]["length"])))
        fig.add_shape(type="rect", x0=extent[0], x1=extent[1],
            y0=lane - .035, y1=lane + .035, line_width=0, fillcolor=p[color])
        for gene in genes:
            if gene["accession"] != accession:
                continue
            start, end = sorted((display(accession, gene["start"]), display(accession, gene["end"])))
            if window and (end <= window[0] or start >= window[1]):
                continue
            y = lane + (.13 if lane else -.13)
            optix = _is_optix(gene)
            fig.add_scatter(x=[start, end], y=[y, y], mode="lines", showlegend=False,
                line={"color": p["red"] if optix else p["muted"], "width": 7 if optix else 4},
                text=[html.escape(str(gene["name"])) + f"<br>[{gene['start']:,}, {gene['end']:,})"] * 2,
                hovertemplate="%{text}<extra></extra>")
            strand = _gene_strand(gene["strand"]) * signs[accession]
            fig.add_scatter(x=[end if strand == 1 else start], y=[y], mode="markers", showlegend=False,
                marker={"symbol": "triangle-right" if strand == 1 else "triangle-left", "size": 10,
                        "color": p["red"] if optix else p["muted"]}, hoverinfo="skip")
            if optix:
                fig.add_annotation(x=(start + end) / 2, y=y + (.11 if lane else -.11),
                    text="<i>optix</i>", showarrow=False, font={"size": 16, "color": p["red"]})
    title = "Read one alignment as a ribbon" if first_only else "Connect the matching DNA intervals"
    subtitle = f"{len(selected):,} measured alignment{'s' if len(selected) != 1 else ''}; transparent bands connect the reported endpoints"
    if centered:
        title = "Place optix at the same display origin"
        subtitle = ("Query reflected to match gene orientation; both optix midpoints are zero" if signs[query] == -1
                    else "Both optix midpoints are zero; original coordinates remain in hover details")
    _layout(fig, title, theme, height=600, subtitle=subtitle, left=185)
    fig.update_layout(margin={"b": 140}, legend={"y": -.29})
    fig.update_yaxes(range=[-.4, 1.4], showgrid=False, tickvals=[0, 1],
        ticktext=[_label(by_id[query]), _label(by_id[target])], fixedrange=True, tickfont={"size": 13})
    fig.update_xaxes(title="Distance from each optix midpoint (bp)" if centered else "Position in imported interval (bp)",
                    range=list(window) if window else None)
    return fig


def optix_gene_models(genes: list[Record], exons: list[Record], hits: list[Record], theme="dark") -> go.Figure:
    """Show the fetched exon intervals, retaining each accession's gene orientation."""
    p = PALETTES[theme]
    target, query = hits[0]["target_accession"], hits[0]["query_accession"]
    centers, signs = _display_transforms(genes, target, query, True)
    fig = go.Figure()
    counts = []
    for accession, lane, color in ((target, 1, "cyan"), (query, 0, "mint")):
        gene = next(g for g in genes if g["accession"] == accession and _is_optix(g))
        selected = [e for e in exons if e["accession"] == accession]
        if not selected:
            raise ValueError(f"No annotated optix exons for {accession}")
        for exon in selected:
            if not gene["start"] <= exon["start"] < exon["end"] <= gene["end"]:
                raise ValueError("An optix exon lies outside its annotated gene interval")
            if _gene_strand(exon["strand"]) != _gene_strand(gene["strand"]):
                raise ValueError("Optix exon and gene strand annotations disagree")
        counts.append(len(selected))
        display = lambda x: signs[accession] * (x - centers[accession])
        fig.add_scatter(x=[display(gene["start"]), display(gene["end"])], y=[lane, lane], mode="lines",
            line={"color": p[color], "width": 2}, name=accession, showlegend=False,
            hovertemplate=f"{html.escape(str(gene['name']))}<br>Gene span: [{gene['start']:,}, {gene['end']:,})<extra></extra>")
        for exon in selected:
            start, end = sorted((display(exon["start"]), display(exon["end"])))
            fig.add_shape(type="rect", x0=start, x1=end, y0=lane-.09, y1=lane+.09,
                line_width=0, fillcolor=p[color])
            fig.add_scatter(x=[(start+end)/2], y=[lane], mode="markers", marker={"size": 15,"opacity":0},
                showlegend=False, text=[f"{html.escape(str(exon['name']))}<br>Exon: [{exon['start']:,}, {exon['end']:,})<br>{exon['end']-exon['start']:,} bp"],
                hovertemplate="%{text}<extra></extra>")
        end = display(gene["end"] if _gene_strand(gene["strand"]) == 1 else gene["start"])
        arrow = "triangle-right" if _gene_strand(gene["strand"]) * signs[accession] == 1 else "triangle-left"
        fig.add_scatter(x=[end], y=[lane], mode="markers", marker={"symbol":arrow,"size":12,"color":p[color]},
                       showlegend=False, hoverinfo="skip")
    _layout(fig,"The current annotations describe different gene spans",theme,height=450,
            subtitle="Boxes are annotated exons; thin lines span the intervening introns",left=180)
    fig.update_yaxes(tickvals=[0,1],ticktext=[f"{query}<br>{counts[1]} annotated exons",f"{target}<br>{counts[0]} annotated exon"],
                     range=[-.5,1.5],showgrid=False,tickfont={"size":13})
    fig.update_xaxes(title="Distance from each optix midpoint (bp)")
    return fig


def locus_context(seq: list[Record], hits: list[Record], windows: list[Record], genes: list[Record],
                  theme="dark", zoom=(-25_000, 125_000)) -> go.Figure:
    """Measured replacement for the old tutorial's unavailable TE/ATAC panels."""
    p = PALETTES[theme]
    target, query = hits[0]["target_accession"], hits[0]["query_accession"]
    shifts, signs = _display_transforms(genes, target, query, True)
    center = ribbons(seq, hits, genes, theme, centered=True, window=zoom)
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, row_heights=[.21, .58, .21],
                        vertical_spacing=.08, subplot_titles=[f"Target {target}: local GC", "Alignment and annotated genes", f"Query {query}: local GC"])
    for accession, row, color in ((target, 1, "cyan"), (query, 3, "mint")):
        def in_zoom(w):
            lo, hi = sorted((signs[accession] * (w["start"] - shifts[accession]), signs[accession] * (w["end"] - shifts[accession])))
            return hi > zoom[0] and lo < zoom[1]
        chosen = sorted((w for w in windows if w["accession"] == accession and in_zoom(w)), key=lambda w:signs[accession] * w["start"])
        if not chosen:
            raise ValueError(f"No composition windows intersect the optix zoom for {accession}")
        fig.add_scatter(x=[signs[accession] * ((w["start"] + w["end"]) / 2 - shifts[accession]) for w in chosen],
            y=[w["gc_percent"] for w in chosen], mode="lines", line={"color": p[color], "width": 2},
            customdata=[[w["start"], w["end"], w["n_fraction"] * 100] for w in chosen], showlegend=False,
            hovertemplate="[%{customdata[0]:,}, %{customdata[1]:,})<br>GC: %{y:.2f}%<br>N: %{customdata[2]:.2f}%<extra></extra>", row=row,col=1)
        fig.update_yaxes(title="GC (%)", range=[0, 100], row=row, col=1)
    for trace in center.data:
        fig.add_trace(trace, row=2, col=1)
    for shape in center.layout.shapes:
        data = shape.to_plotly_json()
        data.update(xref="x2", yref="y2")
        fig.add_shape(data)
    for annotation in center.layout.annotations:
        data = annotation.to_plotly_json()
        data.update(xref="x2", yref="y2")
        fig.add_annotation(data)
    fig.update_yaxes(range=[-.4, 1.4], showgrid=False, tickvals=[0, 1], ticktext=["Query", "Target"], row=2, col=1)
    _layout(fig, "The optix neighborhood: sequence and alignment", theme, height=880,
            subtitle="Fresh NCBI sequence/annotation data; GC is measured, not an ATAC-seq assay")
    fig.update_xaxes(range=list(zoom), title="Distance from each optix midpoint (bp)", row=3, col=1)
    fig.update_annotations(font={"size": 15, "color": p["text"]})
    fig.update_layout(legend={"y": -.1})
    return fig


def make_figures(sequences, matches, gc_windows, genes=(), *, theme="dark", optix_exons=()) -> dict[str, go.Figure]:
    """Create figures only from validated query results; optix views need annotations."""
    if theme not in PALETTES:
        raise ValueError("theme must be dark or light")
    seq, hits, windows, annotations = validate_inputs(sequences, matches, gc_windows, genes)
    exons = records(optix_exons)
    figures = {
        "sequence-lengths": sequence_lengths(seq, theme),
        "base-counts": base_counts(seq, theme),
        "gc-content": gc_content(seq, theme),
        "gc-windows": composition_windows(seq, windows, theme),
        "alignment-dot-plot": dot_plot(hits, theme),
        "first-alignment-ribbon": ribbons(seq, hits, (), theme, first_only=True),
        "alignment-ribbons": ribbons(seq, hits, annotations, theme),
    }
    if all("source_length" in row for row in seq):
        figures["source-record-lengths"] = source_record_lengths(seq, theme)
    if annotations:
        # Fail on missing or ambiguous optix instead of silently drawing a wrong origin.
        figures["optix-centered-ribbons"] = ribbons(seq, hits, annotations, theme, centered=True)
        figures["optix-neighborhood"] = locus_context(seq, hits, windows, annotations, theme)
        if exons:
            figures["optix-gene-models"] = optix_gene_models(annotations, exons, hits, theme)
    return figures


def _interactive_html(dark: go.Figure, light: go.Figure, path: Path) -> None:
    """Write a local Plotly document with automatic theme and an explicit toggle."""
    import plotly.io as pio
    from plotly.utils import PlotlyJSONEncoder
    config = {"responsive": True, "displaylogo": False,
              "toImageButtonOptions": {"format": "svg", "filename": path.stem}}
    themes = json.dumps({"dark": dark.to_plotly_json(), "light": light.to_plotly_json()}, cls=PlotlyJSONEncoder).replace("</", "<\\/")
    post = """
const plots = __THEMES__;
const graph = document.getElementById('{plot_id}');
const button = document.createElement('button');
button.type = 'button'; button.className = 'theme-toggle';
document.body.prepend(button);
const preference = window.matchMedia('(prefers-color-scheme: light)');
let mode = preference.matches ? 'light' : 'dark';
let manual = false;
function apply(theme) {
  mode = theme; document.documentElement.dataset.theme = theme;
  button.textContent = theme === 'light' ? 'Switch to dark mode' : 'Switch to light mode';
  button.setAttribute('aria-pressed', String(theme === 'light'));
  const layout = structuredClone(plots[theme].layout);
  layout.autosize = true; layout.width = null; layout.title = {text: ''};
  layout.margin.t = 70;
  if (graph.clientWidth < 666) {
    layout.font.size = 12;
    layout.margin.l = Math.min(layout.margin.l, 170);
    layout.margin.r = 18;
    layout.legend.font.size = 11;
    for (const key of Object.keys(layout)) {
      if (/^[xy]axis[0-9]*$/.test(key)) {
        layout[key].tickfont = {...layout[key].tickfont, size: 11};
        if (layout[key].title) layout[key].title.font = {...layout[key].title.font, size: 12};
        layout[key].nticks = 5;
      }
    }
    for (const annotation of layout.annotations || []) {
      annotation.font = {...annotation.font, size: 11};
    }
  }
  Plotly.react(graph, plots[theme].data, layout, __CONFIG__);
}
button.addEventListener('click', () => { manual = true; apply(mode === 'light' ? 'dark' : 'light'); });
preference.addEventListener('change', () => { if (!manual) apply(preference.matches ? 'light' : 'dark'); });
let resizeFrame;
window.addEventListener('resize', () => {
  cancelAnimationFrame(resizeFrame);
  resizeFrame = requestAnimationFrame(() => apply(mode));
});
apply(mode);
""".replace("__THEMES__", themes).replace("__CONFIG__", json.dumps(config))
    interactive = go.Figure(dark)
    interactive.update_layout(width=None, autosize=True, title={"text": ""})
    pio.write_html(interactive, str(path), include_plotlyjs="directory", full_html=True,
                   post_script=post, config=config, auto_open=False, default_width="100%")
    source = path.read_text()
    title_parts = dark.layout.title.text.split("<br>", 1)
    heading = "<header><h1>" + title_parts[0] + "</h1>"
    if len(title_parts) == 2:
        heading += "<p>" + title_parts[1].replace("<sup>", "").replace("</sup>", "") + "</p>"
    heading += "</header>"
    style = """<meta name="viewport" content="width=device-width, initial-scale=1">
<style>*{box-sizing:border-box}html{color-scheme:dark;background:#071320;color:#eef4f9}html[data-theme=light]{color-scheme:light;background:#f5f8fc;color:#173247}body{margin:0 auto;max-width:1600px;padding:18px;font-family:Arial,sans-serif}header{padding:12px 20px 0}h1{font-size:26px;margin:8px 0;line-height:1.25}header p{font-size:15px;line-height:1.5;margin:8px 0}.theme-toggle{display:block;margin:0 0 10px auto;padding:9px 14px;background:transparent;color:inherit;border:1px solid currentColor;border-radius:7px;cursor:pointer}.plotly-graph-div{width:100%!important;min-height:420px}@media(max-width:666px){h1{font-size:22px}header{padding:8px 10px 0}}@media(max-width:568px){body{padding:6px}h1{font-size:20px}.plotly-graph-div{min-height:460px}}</style>
"""
    source = source.replace("<head>", "<head>" + style, 1)
    source = source.replace("<body>", "<body>" + heading, 1)
    path.write_text(source)


def write_figures(sequences, matches, gc_windows, genes, output_dir: str | Path, *,
                  svg_transform: Callable[[str], str] | None = None, optix_exons=()) -> dict[str, str]:
    """Export dark/light SVGs, interactive HTML and JSON from the same figures.

    Static export needs Kaleido and Chrome. The website can pass its existing
    theme_svg transformer to create adaptive `<img>` assets. The examples also
    retain explicit light SVGs, so they work independently of that website.
    """
    seq, hits, windows, annotations = validate_inputs(sequences, matches, gc_windows, genes)
    dark = make_figures(seq, hits, windows, annotations, theme="dark", optix_exons=optix_exons)
    light = make_figures(seq, hits, windows, annotations, theme="light", optix_exons=optix_exons)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    # One browser session exports the batch. Starting Chrome once per figure
    # is much slower, and can leave a stalled renderer during a large rerun.
    import plotly.io as pio
    export_figures, export_paths = [], []
    for name in dark:
        export_figures.extend((dark[name], light[name]))
        export_paths.extend((out / f"{name}.svg", out / f"{name}-light.svg"))
    pio.write_images(export_figures, export_paths, format="svg",
                     width=[int(figure.layout.width) for figure in export_figures],
                     height=[int(figure.layout.height) for figure in export_figures])
    paths = {}
    for name, figure in dark.items():
        raw = (out / f"{name}.svg").read_text()
        title = figure.layout.title.text.split("<br>")[0]
        at = raw.index(">")
        raw = raw[:at] + ' role="img" aria-label="' + html.escape(title, quote=True) + '"' + raw[at:]
        (out / f"{name}.svg").write_text(svg_transform(raw) if svg_transform else raw)
        figure.write_json(out / f"{name}.json", pretty=True)
        _interactive_html(figure, light[name], out / f"{name}.html")
        paths[name] = str(out / f"{name}.svg")
    return paths
