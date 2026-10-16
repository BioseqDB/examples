"""Read checksum-pinned NCBI gene models without changing their assembly coordinates."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
from urllib.parse import unquote
from urllib.request import Request, urlopen


def _download_gff(record: dict, cache_dir: Path) -> Path:
    """Cache the annotation, checking its checksum on every use."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{record['assembly']}_genomic.gff.gz"
    if not path.exists():
        request = Request(record["gff_url"], headers={"User-Agent": "BioseqDB-butterfly-tutorial/1.0"})
        with urlopen(request, timeout=60) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != record["gff_sha256"]:
            raise ValueError(f"NCBI annotation checksum changed for {record['assembly']}; review the source before continuing.")
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_bytes(data)
        temporary.replace(path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["gff_sha256"]:
        raise ValueError(f"Cached annotation checksum mismatch: {path}. Remove it and retry the download.")
    return path


def _features(record: dict, cache_dir: Path):
    with gzip.open(_download_gff(record, cache_dir), "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError("Expected nine columns in the pinned NCBI GFF3 file.")
            if fields[0] != record["accession"]:
                continue
            start0, end0 = int(fields[3]) - 1, int(fields[4])
            if start0 >= record["end0"] or end0 <= record["start0"]:
                continue
            attributes = {
                key: unquote(value)
                for key, value in (part.split("=", 1) for part in fields[8].split(";") if "=" in part)
            }
            yield fields[2], start0, end0, fields[6], fields[7], attributes


def _checked_optix_features(record: dict, features: list) -> list:
    """Check the narrow, single-transcript models used in this comparison.

    This tutorial does not implement arbitrary GFF translation exceptions or
    isoform selection. Reject those cases rather than extracting an unrelated
    CDS merely because its coordinates were entered in the manifest.
    """
    gene = record["optix_gene"]
    genes = [feature for feature in features
             if feature[0] == "gene" and feature[5].get("ID") == gene["id"]]
    if len(genes) != 1 or genes[0][1:4] != (gene["start0"], gene["end0"], gene["strand"]):
        raise ValueError(f"The pinned optix model disagrees with sources.json for {record['accession']}.")
    transcripts = [feature for feature in features
                   if feature[0] in {"mRNA", "transcript"}
                   and gene["id"] in feature[5].get("Parent", "").split(",")]
    if len(transcripts) != 1:
        raise ValueError("This optix comparison requires exactly one annotated transcript per gene.")
    transcript = transcripts[0]
    transcript_id = transcript[5]["ID"]
    if transcript[5].get("Parent") != gene["id"]:
        raise ValueError("The optix transcript must belong only to the pinned gene.")
    children = [feature for feature in features
                if feature[0] in {"exon", "CDS"}
                and transcript_id in feature[5].get("Parent", "").split(",")]
    for kind, start0, end0, strand, phase, attrs in genes + transcripts + children:
        if any(key in attrs for key in ("exception", "transl_except", "start_range", "end_range")) or attrs.get("partial") == "true":
            raise ValueError("Partial or exceptional optix models need a separate translation workflow.")
        if not gene["start0"] <= start0 < end0 <= gene["end0"] or strand != gene["strand"]:
            raise ValueError("An optix feature disagrees with the gene interval or strand.")
        if kind in {"exon", "CDS"}:
            if attrs.get("Parent") != transcript_id or not transcript[1] <= start0 < end0 <= transcript[2]:
                raise ValueError("An optix exon or CDS disagrees with its transcript parent.")
        if kind == "CDS" and phase != "0":
            raise ValueError("This optix comparison requires the pinned phase-zero CDS model.")
    for kind, manifest_key in (("exon", "optix_exons"), ("CDS", "optix_cds")):
        actual = sorted((start0, end0, strand, attrs.get("protein_id") if kind == "CDS" else None)
                        for feature_kind, start0, end0, strand, _, attrs in children if feature_kind == kind)
        expected = sorted((item["start0"], item["end0"], item["strand"], item.get("protein_id") if kind == "CDS" else None)
                          for item in record[manifest_key])
        if not actual or actual != expected:
            raise ValueError(f"The pinned optix {kind} features disagree with sources.json for {record['accession']}.")
    return children


def load_annotations(records: list[dict], cache_dir: Path = Path("output/downloads")) -> list[dict]:
    """Return real genes/pseudogenes overlapping each sequence window.

    Plot coordinates ``start`` and ``end`` are 0-based, half-open and relative
    to the downloaded window, clipped at its edges. Absolute coordinates and
    clipping flags remain available for hover labels. An optix label is only
    assigned to the exact GFF ID recorded in sources.json.
    """
    genes = []
    for record in records:
        features = list(_features(record, Path(cache_dir)))
        _checked_optix_features(record, features)
        matched_optix = False
        for kind, start0, end0, strand, _, attrs in features:
            if kind not in {"gene", "pseudogene"}:
                continue
            is_optix = attrs.get("ID") == record["optix_gene"]["id"]
            if is_optix:
                expected = record["optix_gene"]
                if (start0, end0, strand) != (expected["start0"], expected["end0"], expected["strand"]):
                    raise ValueError(f"The pinned optix model disagrees with sources.json for {record['accession']}.")
                matched_optix = True
            genes.append({
                "accession": record["accession"],
                "start": max(start0, record["start0"]) - record["start0"],
                "end": min(end0, record["end0"]) - record["start0"],
                "name": attrs.get("gene") or attrs.get("Name") or attrs["ID"],
                "strand": strand,
                "is_optix": is_optix,
                "feature_type": kind,
                "id": attrs["ID"],
                "description": attrs.get("description", ""),
                "absolute_start0": start0,
                "absolute_end0": end0,
                "clipped": start0 < record["start0"] or end0 > record["end0"],
            })
        if not matched_optix:
            raise ValueError(f"No pinned optix gene model found in the selected {record['accession']} window.")
    return genes


def load_optix_exons(records: list[dict], cache_dir: Path = Path("output/downloads")) -> list[dict]:
    """Return only exons belonging to the pinned gene's transcripts."""
    result = []
    for record in records:
        features = list(_features(record, Path(cache_dir)))
        children = _checked_optix_features(record, features)
        for kind, start0, end0, strand, _, attrs in children:
            if kind != "exon":
                continue
            result.append({
                "accession": record["accession"],
                "start": max(start0, record["start0"]) - record["start0"],
                "end": min(end0, record["end0"]) - record["start0"],
                "name": attrs["ID"],
                "strand": strand,
                "is_optix": True,
            })
    return result


def coding_sequence(window_sequence: str, record: dict) -> str:
    """Extract the manifest's CDS in transcript orientation from the SQL result."""
    if len(window_sequence) != record["end0"] - record["start0"]:
        raise ValueError("Sequence length does not match the pinned window.")
    intervals = sorted(record["optix_cds"], key=lambda item: item["start0"])
    strand = record["optix_gene"]["strand"]
    if strand not in {"+", "-"} or any(item["strand"] != strand for item in intervals):
        raise ValueError("CDS strand metadata is inconsistent.")
    pieces = []
    for item in intervals:
        if not record["start0"] <= item["start0"] < item["end0"] <= record["end0"]:
            raise ValueError("A coding interval falls outside the downloaded window.")
        pieces.append(window_sequence[item["start0"] - record["start0"]:item["end0"] - record["start0"]])
    sequence = "".join(pieces).upper()
    if any(base not in "ACGT" for base in sequence):
        raise ValueError("This coding-sequence comparison expects only unambiguous DNA.")
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1] if strand == "-" else sequence


def compare_coding_sequences(erato: str, melpomene: str) -> dict:
    """Compare equal-length annotated CDSs; this is not a general aligner.

    The two pinned CDSs share a reading frame and have no length difference.
    Reject other inputs rather than silently treating indels as substitutions.
    Protein counts exclude their shared terminal stop codon.
    """
    if len(erato) != len(melpomene) or not erato or len(erato) % 3:
        raise ValueError("Expected nonempty CDSs of the same length, divisible by three.")
    bases = "TCAG"
    amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    code = dict(zip((a + b + c for a in bases for b in bases for c in bases), amino_acids))
    proteins = []
    for sequence in (erato, melpomene):
        if any(base not in "ACGT" for base in sequence):
            raise ValueError("Expected unambiguous uppercase coding sequences.")
        protein = "".join(code[sequence[i:i + 3]] for i in range(0, len(sequence), 3))
        if not sequence.startswith("ATG") or not protein.endswith("*") or "*" in protein[:-1]:
            raise ValueError("Expected a complete CDS with one terminal stop codon.")
        proteins.append(protein[:-1])
    return {
        "cds_length_bp": len(erato),
        "identical_nucleotides": sum(a == b for a, b in zip(erato, melpomene)),
        "protein_length_amino_acids": len(proteins[0]),
        "identical_amino_acids": sum(a == b for a, b in zip(*proteins)),
        "amino_acid_differences": [
            {"position1": i + 1, "erato": a, "melpomene": b}
            for i, (a, b) in enumerate(zip(*proteins)) if a != b
        ],
    }
