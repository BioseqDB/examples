#!/usr/bin/env python3
"""Plot the CSV files exported by sql/09_export.sql.

All alignment, translation, counting and normalization happen in PostgreSQL.
This script reads those results and draws transparent SVGs for both site themes.
Run: python render_charts.py --results results --output figures
"""
from __future__ import annotations

import argparse
import csv
from html import escape
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, to_hex
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter, MaxNLocator, FuncFormatter
import numpy as np

ROOT = Path(__file__).resolve().parent
AA = list("ARNDCQEGHILKMFPSTWYV")
DNA = list("ACGT")
COLORS = {"A": "#74c8a8", "C": "#65b9dd", "G": "#f4894f", "T": "#c9a8ff"}
TEXT, MUTED, GRID = "#eef4f9", "#a9bccb", "#28425b"
LIGHT = {TEXT: "#173247", MUTED: "#405c72", GRID: "#bfcedc",
         "#74c8a8": "#167251", "#65b9dd": "#146c91",
         "#f4894f": "#b44d1c", "#c9a8ff": "#735091",
         "#e0b04a": "#805b16", "#20364a": "#ccd9e4", "#071320": "#f5f8fc"}
SCORE_MAP = LinearSegmentedColormap.from_list("bioseq_scores", ["#f4894f", "#20364a", "#65b9dd"], N=65)
RATIO_MAP = LinearSegmentedColormap.from_list("bioseq_ratios", ["#20364a", "#65b9dd", "#74c8a8"], N=65)
# Map all sampled colors, including those used in the colorbars. Mapping only
# the three anchors leaves dark intermediate colors and creates false bands.
LIGHT_SCORE_MAP = LinearSegmentedColormap.from_list(
    "bioseq_scores_light", ["#b44d1c", "#e8eff5", "#146c91"], N=65)
LIGHT_RATIO_MAP = LinearSegmentedColormap.from_list(
    "bioseq_ratios_light", ["#e8eff5", "#146c91", "#167251"], N=65)
for dark_map, light_map in ((SCORE_MAP, LIGHT_SCORE_MAP), (RATIO_MAP, LIGHT_RATIO_MAP)):
    for sample in np.linspace(0, 1, dark_map.N):
        LIGHT[to_hex(dark_map(sample))] = to_hex(light_map(sample))
SCORE_MAP.set_bad("#28425b")
RATIO_MAP.set_bad("#28425b")
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 13, "text.color": TEXT,
    "axes.labelcolor": TEXT, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.edgecolor": GRID, "grid.color": GRID, "axes.spines.top": False,
    "axes.spines.right": False, "axes.spines.left": False,
    "axes.spines.bottom": False, "svg.fonttype": "none",
    "svg.hashsalt": "bioseq-scoring-matrices", "savefig.transparent": True,
})


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def number(value):
    return float(value) if value not in (None, "", "NULL") else np.nan


def matrix(records, alphabet, x, y, value, divisor=1):
    result = np.full((len(alphabet), len(alphabet)), np.nan)
    indices = {letter: i for i, letter in enumerate(alphabet)}
    for row in records:
        result[indices[row[x]], indices[row[y]]] = number(row[value]) / divisor
    return result


def figure(width=12, height=7):
    fig, ax = plt.subplots(figsize=(width, height), dpi=72)
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    ax.tick_params(length=0, pad=9)
    return fig, ax


def save(fig, path, title, description):
    fig.savefig(path, format="svg", metadata={"Date": None, "Creator": "BioseqDB scoring_matrices/render_charts.py"})
    plt.close(fig)
    source = path.read_text().replace("'DejaVu Sans'", "Arial, sans-serif")
    source = re.sub(r'width="([0-9.]+)pt" height="([0-9.]+)pt"', r'width="\1" height="\2"', source, count=1)
    for dark in LIGHT:
        for prop in ("fill", "stroke"):
            source = source.replace(f"{prop}: {dark}", f"{prop}: var(--theme-{dark[1:]}, {dark})")
    used = sorted(set(re.findall(r"var\((--theme-[0-9a-f]+),", source)))
    declarations = "".join(f"{token}:{LIGHT['#' + token.removeprefix('--theme-')]};" for token in used)
    metadata = ('\n<style id="bioseq-light-palette">@media(prefers-color-scheme:light){:root{'
                + declarations + '}}</style>\n<title id="chart-title">' + escape(title)
                + '</title>\n<desc id="chart-description">' + escape(description) + '</desc>\n')
    source = source.replace("<svg ", '<svg data-bioseq-palette="explicit" role="img" aria-labelledby="chart-title chart-description" ', 1)
    start = source.index(">", source.index("<svg ")) + 1
    source = source[:start] + metadata + source[start:]
    path.write_text(source)


def heatmap(values, alphabet, output, title, subtitle, *, ratios=False, annotate=False,
            xlabel="Partner amino acid", ylabel="First amino acid", small=False):
    fig, ax = figure(11, 10.5 if not small else 8.8)
    fig.text(.09, .96, title, fontsize=22, weight="bold", va="top")
    fig.text(.09, .915, subtitle, fontsize=12, color=MUTED, va="top")
    fig.subplots_adjust(left=.13, right=.85, bottom=.13, top=.85)
    kwargs = {"cmap": RATIO_MAP, "vmin": 0, "vmax": 1} if ratios else {
        "cmap": SCORE_MAP, "norm": TwoSlopeNorm(vmin=-24, vcenter=0, vmax=8)}
    edges = np.arange(len(alphabet)+1)-.5
    im = ax.pcolormesh(edges, edges, np.ma.masked_invalid(values), shading="flat", **kwargs)
    ax.set_xlim(-.5, len(alphabet)-.5)
    ax.set_ylim(len(alphabet)-.5, -.5)
    ax.set_aspect("equal")
    ax.set_xticks(range(len(alphabet)), alphabet, fontsize=17)
    ax.set_yticks(range(len(alphabet)), alphabet, fontsize=17)
    ax.set_xlabel(xlabel, labelpad=14)
    ax.set_ylabel(ylabel, labelpad=14)
    ax.set_xticks(np.arange(-.5, len(alphabet)), minor=True)
    ax.set_yticks(np.arange(-.5, len(alphabet)), minor=True)
    ax.grid(which="minor", color=GRID, linewidth=.6)
    ax.tick_params(which="minor", length=0)
    for i in range(len(alphabet)):
        for j in range(len(alphabet)):
            value = values[i, j]
            if np.isnan(value):
                ax.plot([j-.28, j+.28], [i+.28, i-.28], color=MUTED, lw=1.2, alpha=.8)
            elif annotate:
                label = f"{100*value:.1f}%" if ratios else f"{value:.2f}"
                ax.text(j, i, label, ha="center", va="center", fontsize=19,
                        bbox={"boxstyle": "round,pad=.3", "facecolor": "#071320",
                              "edgecolor": "none", "alpha": .86})
    cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.04, shrink=.86)
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor("face")
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0)
    if ratios:
        cb.ax.yaxis.set_major_formatter(PercentFormatter(1))
        cb.set_label("Fraction among differing pairs", labelpad=10)
    else:
        cb.set_label("Log odds (bits)", labelpad=10)
    fig.text(.13, .035, "Slash = excluded diagonal or undefined row" if ratios else
             "Slash = unobserved pair; its unsmoothed score is −∞", color=MUTED, fontsize=11)
    save(fig, output, title, subtitle + ". Axes use the same ordered alphabet. "
         + ("Each off-diagonal row sums to one." if ratios else "Positive scores indicate enrichment relative to background frequencies."))


def render(results, references, output):
    output.mkdir(parents=True, exist_ok=True)
    nucleotide = rows(results / "nucleotide_matrix.csv")
    protein = rows(results / "protein_matrix.csv")
    spectrum = rows(results / "reference_spectrum.csv")
    for name, label, scale in (("blosum62", "BLOSUM62", 2), ("pam250", "PAM250", 3)):
        records = rows(references / f"{name}_scores.csv")
        values = matrix(records, AA, "aa_from", "aa_to", "score", scale)
        heatmap(values, AA, output / f"{name}.svg", label,
                f"NCBI score table · {scale} integer units per bit · same color scale in both reference plots")
    heatmap(matrix(nucleotide, DNA, "base_from", "base_to", "score_half_bits", 2), DNA,
            output / "nucleotide-scores.svg", "DNA log odds from the cohort",
            "Aligned A/C/G/T columns · all distinct sample pairs · N and gaps excluded",
            annotate=True, small=True, xlabel="Partner nucleotide", ylabel="First nucleotide")
    heatmap(matrix(nucleotide, DNA, "base_from", "base_to", "fraction_of_differences"), DNA,
            output / "nucleotide-partners.svg", "Partner bases across all sample pairs",
            "Each row contains differing pairs only · all 1,000 samples contribute",
            ratios=True, annotate=True, small=True,
            xlabel="Partner nucleotide", ylabel="First nucleotide")
    heatmap(matrix(protein, AA, "aa_from", "aa_to", "score_half_bits", 2), AA,
            output / "protein-scores.svg", "Protein log odds from the cohort",
            "13 mitochondrial coding regions · accepted codons · same bit scale as the reference plots")
    heatmap(matrix(protein, AA, "aa_from", "aa_to", "fraction_of_differences"), AA,
            output / "protein-partners-cohort.svg", "Which amino acids differ in the cohort?",
            "Partner frequencies among differing pairs · defined rows sum to 100%", ratios=True)
    baseline = [r for r in rows(results / "reference_protein_matrix.csv") if r["matrix"] == "BLOSUM62"]
    heatmap(matrix(baseline, AA, "aa_from", "aa_to", "fraction_of_differences"), AA,
            output / "protein-partners-blosum62.svg", "The same comparison in BLOSUM62",
            "Original NCBI joint probabilities · identical denominator and color scale", ratios=True)

    fig, ax = figure(11, 7)
    fig.subplots_adjust(left=.13, right=.95, top=.8, bottom=.16)
    fig.text(.08, .96, "Where reference-base differences go", fontsize=22, weight="bold", va="top")
    fig.text(.08, .89, "Each bar includes differences only; unchanged bases are outside this denominator.", fontsize=12, color=MUTED)
    for i, origin in enumerate(DNA):
        left = 0
        for target in DNA:
            if target == origin:
                continue
            row = next(r for r in spectrum if r["reference_base"] == origin and r["sample_base"] == target)
            fraction = number(row["fraction_of_differences"])
            ax.barh(i, fraction, left=left, height=.56, color=COLORS[target])
            if fraction > .09:
                ax.text(left + fraction/2, i, f"{target}\n{fraction:.1%}", va="center", ha="center", fontsize=15, color="#071320")
            left += fraction
    ax.set_yticks(range(4), [f"{a} →" for a in DNA], fontsize=19)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(PercentFormatter(1))
    ax.set_xlabel("Fraction of differences from this reference base")
    ax.legend(handles=[Patch(color=COLORS[b], label=b) for b in DNA], ncol=4,
              frameon=False, loc="upper center", bbox_to_anchor=(.5, -.16))
    save(fig, output / "nucleotide-destinations.svg", "Reference-to-sample nucleotide differences",
         "Four stacked bars show the fraction of each alternate base among callable differences. Reference direction does not identify the ancestral mutation direction.")

    positions = rows(results / "position_summary.csv")
    # Coordinates and fractions have already been counted in SQL. Display raw
    # per-position counts without smoothing away isolated differences.
    x = [int(r["position1"]) for r in positions]
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), dpi=72, sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})
    fig.subplots_adjust(left=.1, right=.97, top=.83, bottom=.13, hspace=.26)
    fig.text(.07, .96, "Differences and missing calls along the mitochondrial genome", fontsize=20, weight="bold", va="top")
    fig.text(.07, .90, "Counts at each rCRS position · individual sites, with no window averaging", color=MUTED, fontsize=12)
    changed = [r for r in positions if int(r["reference_differences"]) > 0]
    axes[0].vlines([int(r["position1"]) for r in changed], 0,
                  [int(r["reference_differences"]) for r in changed], color=COLORS["A"], lw=.8)
    axes[0].set_ylabel("Samples differing\nfrom the reference")
    axes[1].plot(x, [int(r["unknown_n"]) for r in positions], color=COLORS["G"], lw=.8, label="Sample N")
    axes[1].plot(x, [int(r["gaps"]) for r in positions], color=COLORS["C"], lw=.8, label="Deletion / uncovered position")
    axes[1].set_ylabel("Excluded calls")
    axes[1].set_xlabel("Reference coordinate (1-based)")
    axes[1].legend(frameon=False, loc="lower right", bbox_to_anchor=(1, 1.04), fontsize=10)
    axes[1].xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
    for ax in axes:
        ax.set_axisbelow(True)
        ax.grid(axis="y", alpha=.5)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    save(fig, output / "genome-differences.svg", "Differences and missing calls along the mitochondrial genome",
         "The upper panel counts sample bases differing from canonical rCRS bases at each position. The lower panel separately counts N and deletions or uncovered positions.")

    quality = [r for r in rows(results / "genome_quality.csv") if r["role"] == "sample"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), dpi=72)
    fig.subplots_adjust(left=.08, right=.96, top=.76, bottom=.2, wspace=.3)
    fig.text(.06, .96, "Check the records before counting substitutions", fontsize=22, weight="bold", va="top")
    fig.text(.06, .88, f"{len(quality):,} complete mitochondrial records · one point per record", fontsize=12, color=MUTED)
    indices = np.arange(1, len(quality)+1)
    axes[0].scatter(indices, [int(r["bp"]) for r in quality], s=8, alpha=.6, color=COLORS["C"], rasterized=False)
    axes[0].set_ylabel("Sequence length (bp)")
    axes[1].scatter(indices, [int(r["n_count"]) for r in quality], s=8, alpha=.6, color=COLORS["G"], rasterized=False)
    axes[1].set_ylabel("Uncalled nucleotides (N)")
    for ax in axes:
        ax.set_xlabel("Record, sorted by accession")
        ax.grid(axis="y", alpha=.35)
        ax.set_axisbelow(True)
    save(fig, output / "record-quality.svg", "Length and N content of the mitochondrial records",
         f"Sequence length and raw N count are shown separately for all {len(quality):,} sample records. The reference is excluded from these panels.")

    cds = rows(results / "cds_quality.csv")
    fig, ax = figure(12, 7)
    fig.subplots_adjust(left=.12, right=.96, top=.76, bottom=.18)
    fig.text(.07, .96, "Coding regions retained for protein comparisons", fontsize=21, weight="bold", va="top")
    fig.text(.07, .89, "Exclude a CDS for gaps, missing alignment, a premature stop or an unusable start.", color=MUTED, fontsize=12)
    x = np.arange(len(cds))
    retained = [int(r["retained_cds"]) for r in cds]
    excluded = [int(r["sample_cds"]) - int(r["retained_cds"]) for r in cds]
    ax.bar(x, retained, color=COLORS["A"], label="Retained")
    ax.bar(x, excluded, bottom=retained, color=COLORS["G"], label="Excluded CDS (count above)")
    # A handful of exclusions is too small to distinguish on a 0–1000 axis.
    # Label the count without truncating the scale or exaggerating the loss.
    for index, count in enumerate(excluded):
        if count:
            ax.annotate(str(count), (index, retained[index] + count),
                        xytext=(0, 5), textcoords="offset points", ha="center",
                        va="bottom", fontsize=12, color=COLORS["G"], weight="bold")
    ax.set_xticks(x, [r["gene"] for r in cds], rotation=45, ha="right", fontsize=12)
    ax.set_ylabel("Sample coding regions")
    ax.set_axisbelow(True)
    ax.grid(axis="y", alpha=.35)
    ax.legend(frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(.5, 1.10))
    save(fig, output / "protein-coverage.svg", "Coding regions retained for amino-acid counting",
             "Each bar starts with the full cohort. A coding region is excluded for insertions, deletions, incomplete alignment, premature stops, or invalid or ambiguous initiation codons. Ambiguous internal codons are counted separately.")

    distances_path = results / "sample_distances.csv"
    if distances_path.exists():
        records = rows(distances_path)
        accessions = sorted({r["sample_from"] for r in records} | {r["sample_to"] for r in records})
        values = matrix(records, accessions, "sample_from", "sample_to", "p_distance")
        fig, ax = figure(11, 10)
        fig.subplots_adjust(left=.19, right=.86, top=.84, bottom=.24)
        fig.text(.08, .96, "Eight samples, compared at callable positions", fontsize=21, weight="bold", va="top")
        fig.text(.08, .90, "p-distance = differing bases / paired A/C/G/T sites", fontsize=12, color=MUTED)
        im = ax.imshow(values, cmap=RATIO_MAP, vmin=0)
        ax.set_xticks(range(len(accessions)), accessions, rotation=45, ha="right", fontsize=11)
        ax.set_yticks(range(len(accessions)), accessions, fontsize=11)
        for i in range(len(accessions)):
            for j in range(len(accessions)):
                ax.text(j, i, f"{100*values[i,j]:.2f}", ha="center", va="center", fontsize=11,
                        bbox={"boxstyle": "round,pad=.2", "facecolor": "#071320",
                              "edgecolor": "none", "alpha": .86})
        cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.04)
        cb.outline.set_visible(False)
        cb.ax.yaxis.set_major_formatter(PercentFormatter(1))
        cb.set_label("Differing fraction")
        save(fig, output / "sample-distances.svg", "Callable-site p-distances for eight cohort samples",
             "Each cell is a pairwise mismatch percentage, with Ns and gaps removed from that pair's denominator. It is a descriptive distance, not a divergence time.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--references", type=Path, default=ROOT / "reference_matrices")
    parser.add_argument("--output", type=Path, default=ROOT / "figures")
    args = parser.parse_args()
    render(args.results, args.references, args.output)
    print(f"Wrote SVG plots to {args.output}")


if __name__ == "__main__":
    main()
