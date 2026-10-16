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
            yield fields[2], start0, end0, fields[6], attributes


def load_annotations(records: list[dict], cache_dir: Path = Path("output/downloads")) -> list[dict]:
    """Return real genes/pseudogenes overlapping each sequence window.

    Plot coordinates ``start`` and ``end`` are 0-based, half-open and relative
    to the downloaded window, clipped at its edges. Absolute coordinates and
    clipping flags remain available for hover labels. An optix label is only
    assigned to the exact GFF ID recorded in sources.json.
    """
    genes = []
    for record in records:
        matched_optix = False
        for kind, start0, end0, strand, attrs in _features(record, Path(cache_dir)):
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
        transcripts = {
            attrs["ID"] for kind, _, _, _, attrs in features
            if kind in {"mRNA", "transcript"} and record["optix_gene"]["id"] in attrs.get("Parent", "").split(",")
        }
        for kind, start0, end0, strand, attrs in features:
            if kind != "exon" or not (transcripts & set(attrs.get("Parent", "").split(","))):
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
