"""Small explanatory SVGs for the butterfly tutorial.

These drawings are conceptual illustrations, not measurements or reconstructed
specimens. Empirical charts are generated separately by plots.py from SQL rows.
"""
from __future__ import annotations

import argparse
import html
import re
from pathlib import Path

from .plots import PALETTES

P = PALETTES["dark"]


def adaptive_svg(source: str) -> str:
    """Give standalone illustrations paired colors without a website dependency."""
    declarations = []
    for name, dark in P.items():
        key = "--theme-" + dark[1:]
        source = source.replace(dark, f"var({key}, {dark})")
        declarations.append(f"{key}:{PALETTES['light'][name]};")
    css = '<style id="bioseq-light-palette">@media(prefers-color-scheme:light){:root{' + "".join(declarations) + '}}</style>'
    return re.sub(r"(<svg\b[^>]*>)", lambda m: m[1] + css, source, count=1)


class Drawing:
    def __init__(self, title: str, description: str, height=520):
        self.height = height
        self.parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 860 {height}" role="img" aria-labelledby="title desc">',
            f'<title id="title">{html.escape(title)}</title><desc id="desc">{html.escape(description)}</desc>',
            '<defs><marker id="arrow" markerUnits="userSpaceOnUse" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto-start-reverse"><path d="M0,0 L11,6 L0,12" fill="' + P["muted"] + '"/></marker></defs>',
        ]
        self.text(30, 43, title, size=27, weight="bold")

    def text(self, x, y, value, *, size=19, color="text", anchor="start", weight="normal", italic=False):
        lines = str(value).split("\n")
        attr = ' font-style="italic"' if italic else ""
        self.parts.append(f'<text x="{x}" y="{y}" fill="{P.get(color,color)}" font-family="Arial, sans-serif" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}"{attr}>')
        for i, line in enumerate(lines):
            self.parts.append(f'<tspan x="{x}" dy="{0 if i == 0 else size * 1.4}">{html.escape(line)}</tspan>')
        self.parts.append('</text>')

    def rect(self, x, y, w, h, color="panel", *, stroke=None, radius=12, opacity=1):
        self.parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" fill="{P.get(color,color)}" opacity="{opacity}"' +
            (f' stroke="{P.get(stroke,stroke)}" stroke-width="1.5"' if stroke else "") + '/>')

    def line(self, x0, y0, x1, y1, color="muted", *, width=2, arrow=False, dash=False):
        self.parts.append(f'<path d="M{x0},{y0} L{x1},{y1}" fill="none" stroke="{P.get(color,color)}" stroke-width="{width}"' +
            (' marker-end="url(#arrow)"' if arrow else "") + (' stroke-dasharray="7 5"' if dash else "") + '/>')

    def path(self, d, color, *, opacity=1, stroke=None, width=2):
        self.parts.append(f'<path d="{d}" fill="{P.get(color,color)}" opacity="{opacity}"' +
            (f' stroke="{P.get(stroke,stroke)}" stroke-width="{width}"' if stroke else "") + '/>')

    def circle(self, x, y, radius, color="cyan", *, opacity=1):
        self.parts.append(f'<circle cx="{x}" cy="{y}" r="{radius}" fill="{P.get(color,color)}" opacity="{opacity}"/>')

    def note(self, text):
        self.text(30, self.height - 25, text, size=16, color="muted")

    def finish(self):
        return adaptive_svg("\n".join(self.parts) + "</svg>\n")


def dna_alignment():
    d = Drawing("DNA letters become alignment columns", "A, C, G and T denote DNA bases. N denotes an unknown base. A toy alignment shows a match, insertion, deletion and mismatch. This is not a butterfly result.", 560)
    colors = {"A":"purple", "C":"cyan", "G":"orange", "T":"mint", "N":"muted"}
    names = {"A":"adenine", "C":"cytosine", "G":"guanine", "T":"thymine", "N":"unknown"}
    for i, (letter, color) in enumerate(colors.items()):
        x = 52 + i * 163
        d.rect(x, 83, 68, 65, color, opacity=.15)
        d.text(x+34, 128, letter, size=37, anchor="middle", color=color, weight="bold")
        d.text(x+34, 178, names[letter], size=17, anchor="middle")
    d.line(30, 207, 830, 207, "grid")
    d.text(30, 245, "One possible alignment", size=21, weight="bold")
    d.text(30, 295, "Target", size=18)
    d.text(30, 366, "Query", size=18)
    target, query = "AC-GTAAC", "ACTGT-AT"
    for row, letters in enumerate((target, query)):
        y = 264 + 72 * row
        for i, letter in enumerate(letters):
            x = 160 + i * 71
            if letter != "-":
                d.rect(x, y, 52, 46, colors[letter], opacity=.14, radius=6)
            d.text(x+26, y+33, letter, size=29, color=colors.get(letter,"muted"), anchor="middle", weight="bold")
    for x, label, color in ((186,"match","mint"),(328,"insertion","cyan"),(541,"deletion","orange"),(683,"mismatch","red")):
        d.line(x, 389, x, 414, color)
        d.text(x, 440, label, size=17, color=color, anchor="middle")
    d.text(30, 486, "The gap symbol “–” is added by the alignment; it is not a DNA base.", size=18)
    d.note("Illustrative letters only. Insertion/deletion are described relative to the target.")
    return d.finish()


def workflow():
    d = Drawing("The analysis stays connected to its data", "NCBI supplies sequence and annotation records. BioseqDB stores sequence and metadata and aligns the intervals. SQLAlchemy executes SQL; Plotly turns returned summaries and matches into charts.", 535)
    boxes = [(40,100,340,130,"cyan","NCBI","Accession + assembly\nSequence + annotation"),
             (480,100,340,130,"mint","PostgreSQL + BioseqDB","Import and retain provenance\nAlign with SQL"),
             (40,315,340,130,"purple","SQLAlchemy","Execute parameterized queries\nReceive measured result rows"),
             (480,315,340,130,"orange","Plotly","Composition and alignment charts\nInteractive HTML + SVG")]
    for x,y,w,h,c,title,body in boxes:
        d.rect(x,y,w,h,stroke=c)
        d.text(x+22,y+38,title,size=24,color=c,weight="bold")
        d.text(x+22,y+75,body,size=18)
    d.line(392,165,463,165,arrow=True)
    d.text(428,140,"SQL",size=17,anchor="middle",color="muted")
    d.line(650,244,650,282,arrow=True)
    d.line(650,282,210,282)
    d.line(210,282,210,303,arrow=True)
    d.text(430,271,"Search and summaries run in the database",size=17,anchor="middle",color="muted")
    d.line(392,380,463,380,arrow=True)
    d.note("Python sends queries and draws their results. No external aligner is invoked.")
    return d.finish()


def coordinate_systems():
    d = Drawing("Keep one coordinate convention", "An illustrative annotation covering positions three through six in one-based inclusive notation becomes the zero-based half-open interval [2, 6). An imported subinterval has its own origin that must be retained.", 490)
    d.text(30,92,"Same four bases, two ways to write their interval",size=20,color="muted")
    letters="ACGTTACG"
    for i, letter in enumerate(letters):
        x=170+i*72
        d.rect(x,172,58,58,"cyan",opacity=.25 if 2<=i<6 else .07,radius=7)
        d.text(x+29,211,letter,size=27,anchor="middle",color="cyan" if 2<=i<6 else "muted",weight="bold")
        d.text(x+29,149,str(i+1),size=18,anchor="middle")
        d.text(x+29,263,str(i),size=18,anchor="middle")
    d.text(30,149,"One-based",size=18)
    d.text(30,263,"Zero-based",size=18)
    d.line(314,285,590,285,"cyan",width=3)
    d.text(450,323,"[3, 6] inclusive → [2, 6) half-open",size=23,anchor="middle",color="cyan",weight="bold")
    d.rect(30,355,800,78,stroke="grid")
    d.text(50,385,"chromosome position = imported interval start + local position",size=20)
    d.text(50,415,"Length = end − start. The base at end is not included.",size=18,color="muted")
    d.note("The numbers above are a coordinate example, not positions of the optix gene.")
    return d.finish()


def rectangle_coordinates():
    d=Drawing("A rectangle represents a DNA interval", "Plotly rectangles use x0, x1, y0 and y1. Genomic x coordinates represent an interval. Vertical positions are arbitrary drawing lanes, not another genomic measurement.", 470)
    d.line(115,345,774,345,arrow=True)
    d.line(115,345,115,105,arrow=True)
    d.rect(250,145,350,155,"cyan",opacity=.18,radius=0)
    d.circle(250,300,6,"cyan");d.circle(600,145,6,"cyan")
    d.line(250,300,250,345,"muted",dash=True);d.line(600,145,600,345,"muted",dash=True)
    d.text(233,327,"(x0, y0)",size=24,anchor="end",color="cyan")
    d.text(615,139,"(x1, y1)",size=24,color="cyan")
    d.text(425,225,"one sequence interval",size=23,anchor="middle")
    d.text(425,264,"x0 = start       x1 = end",size=20,anchor="middle",color="muted")
    d.text(445,389,"Position in the imported sequence (bp)",size=20,anchor="middle")
    d.note("The y values place the track on the page; they have no biological unit.")
    return d.finish()


def ribbon_coordinates():
    d=Drawing("Four endpoints make one alignment ribbon", "A forward alignment ribbon joins target start/end on an upper lane to query end/start on a lower lane. Reverse alignments connect the opposite query endpoint order. These are conceptual endpoints.", 560)
    d.rect(65,130,730,16,"cyan",radius=3)
    d.rect(65,366,730,16,"mint",radius=3)
    d.path("M235,146 L510,146 L655,366 L330,366 Z","cyan",opacity=.2)
    for x,y,c in ((235,146,"cyan"),(510,146,"cyan"),(330,366,"mint"),(655,366,"mint")):
        d.circle(x,y,6,c)
    d.text(65,106,"TARGET",size=18,color="cyan",weight="bold")
    d.text(235,184,"target_start",size=18,anchor="end",color="cyan")
    d.text(524,184,"target_end",size=18,color="cyan")
    d.text(65,416,"QUERY",size=18,color="mint",weight="bold")
    d.text(330,353,"query_start",size=18,anchor="end",color="mint")
    d.text(671,353,"query_end",size=18,color="mint")
    d.line(310,208,492,208,arrow=True)
    d.line(548,239,598,305,arrow=True)
    d.line(562,324,399,324,arrow=True)
    d.text(430,473,"Forward: [target_start, target_end, query_end, query_start]",size=18,anchor="middle")
    d.note("For reverse alignments, swap the two query endpoints when drawing the polygon.")
    return d.finish()


def alignment_display():
    d=Drawing("Changing the view does not change the alignment", "An optix-centered view subtracts each gene midpoint. If gene strands differ, the query axis is reflected. This aligns gene orientation in the drawing without rewriting source coordinates or claiming a coordinate lift-over.", 545)
    d.text(30,91,"Raw assembly directions",size=20,weight="bold")
    d.line(80,150,760,150,"cyan",width=8)
    d.line(80,228,760,228,"mint",width=8)
    d.line(300,150,390,150,"red",width=12,arrow=True)
    d.line(600,228,510,228,"red",width=12,arrow=True)
    d.text(345,130,"optix →",size=18,anchor="middle",color="red")
    d.text(555,270,"← optix",size=18,anchor="middle",color="red")
    d.line(390,162,510,213,"muted",dash=True)
    d.line(300,162,600,213,"muted",dash=True)
    d.text(30,325,"Reflect the query, then place both gene midpoints at zero",size=20,weight="bold")
    d.line(90,372,770,372,"cyan",width=8)
    d.line(90,432,770,432,"mint",width=8)
    d.line(390,372,470,372,"red",width=12,arrow=True)
    d.line(390,432,470,432,"red",width=12,arrow=True)
    d.line(430,345,430,458,"muted",dash=True)
    d.text(430,485,"0",size=21,anchor="middle")
    d.note("Display only: x_target = t − center_target; x_query = orientation × (q − center_query).")
    return d.finish()


def atac_concept():
    d=Drawing("Where an ATAC-seq track would come from", "Conceptual ATAC-seq illustration. DNA wraps around nucleosomes; accessible regions can be tagged by Tn5. Sequencing and mapping assay fragments generates a measured signal track. No ATAC-seq measurements are used in this tutorial.", 560)
    d.text(30,90,"DNA accessibility requires a separate experiment",size=20,color="muted")
    d.path("M60,195 C130,130 180,255 240,195 S340,135 390,195 S490,255 550,195 S670,135 790,195","none",stroke="cyan",width=5)
    for x in (150,285,670,755):
        d.circle(x,195,31,"muted",opacity=.3)
        d.circle(x-8,187,17,"muted",opacity=.55)
    d.text(210,271,"DNA wrapped around\nnucleosomes",size=18,anchor="middle",color="muted")
    d.text(487,156,"accessible region",size=19,anchor="middle",color="cyan")
    for x in (440,485,530):
        d.path(f"M{x-10},205 L{x},180 L{x+10},205 Z","orange")
    d.text(485,264,"Tn5 tags accessible DNA",size=18,anchor="middle",color="orange")
    d.line(485,285,485,330,arrow=True)
    d.rect(55,352,300,91,stroke="grid")
    d.text(205,385,"Sequence assay fragments",size=19,anchor="middle")
    d.text(205,418,"Map them to the same assembly",size=17,anchor="middle",color="muted")
    d.line(370,396,416,396,arrow=True)
    d.line(435,441,792,441,"muted")
    d.path("M435,438 L465,438 L478,422 L486,436 L515,435 L532,378 L545,435 L592,437 L610,409 L618,437 L666,436 L676,393 L690,436 L792,438 L792,441 L435,441 Z","purple",opacity=.35)
    d.text(615,474,"Illustrative signal only",size=17,anchor="middle",color="purple")
    d.note("This tutorial plots sequence-derived GC, not ATAC-seq, expression or enhancer activity.")
    return d.finish()


def butterflies_optix():
    d=Drawing("Two butterflies, one genomic neighborhood", "Conceptual butterfly silhouettes with red wing marks, labelled Heliconius erato demophoon and Heliconius melpomene rosina. A sequence interval containing optix appears below each. These are illustrations, not specimen photographs or measured wing patterns.", 535)
    for center,color,species in ((230,"mint","H. erato demophoon"),(630,"cyan","H. melpomene rosina")):
        d.path(f"M{center},215 C{center-25},140 {center-145},100 {center-150},165 C{center-155},218 {center-95},242 {center-48},236 C{center-125},247 {center-125},322 {center-66},317 C{center-23},310 {center-5},257 {center},225 Z",color,opacity=.2,stroke=color)
        d.path(f"M{center},215 C{center+25},140 {center+145},100 {center+150},165 C{center+155},218 {center+95},242 {center+48},236 C{center+125},247 {center+125},322 {center+66},317 C{center+23},310 {center+5},257 {center},225 Z",color,opacity=.2,stroke=color)
        d.path(f"M{center-48},191 L{center-118},166 L{center-120},188 L{center-57},215 Z","red",opacity=.85)
        d.path(f"M{center+48},191 L{center+118},166 L{center+120},188 L{center+57},215 Z","red",opacity=.85)
        d.line(center,171,center,278,"text",width=8)
        d.line(center,176,center-20,142,"muted")
        d.line(center,176,center+20,142,"muted")
        d.text(center,357,species,size=21,anchor="middle",italic=True)
        d.line(center-150,413,center+150,413,color,width=9)
        d.line(center-25,413,center+25,413,"red",width=13)
        d.text(center,453,"optix + surrounding DNA",size=18,anchor="middle",color="muted")
    d.note("Conceptual illustration. The figures below compare DNA, not photographed wing patterns.")
    return d.finish()


DIAGRAMS = {
    "butterflies-optix": butterflies_optix,
    "dna-alignment": dna_alignment,
    "sql-workflow": workflow,
    "coordinate-systems": coordinate_systems,
    "rectangle-coordinates": rectangle_coordinates,
    "ribbon-coordinates": ribbon_coordinates,
    "alignment-display-transform": alignment_display,
    "atac-concept": atac_concept,
}


def write_diagrams(output_dir: str | Path) -> dict[str, str]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    result = {}
    for name, create in DIAGRAMS.items():
        path = out / f"{name}.svg"
        path.write_text(create())
        result[name] = str(path)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("output/figures"))
    args = parser.parse_args()
    write_diagrams(args.output)
