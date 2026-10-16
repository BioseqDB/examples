"""Retrieve sequences with SQL, measure them, align them, and render the figures.

No aligner executable is called. Python replays the returned CIGAR strings as an
independent check; BioseqDB performs the search and base-level alignment.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import time

from sqlalchemy import create_engine, text

ROOT = Path(__file__).resolve().parents[1]


def sql(name):
    return text((ROOT / "sql" / name).read_text())


def rows(connection, name, parameters=None):
    return [dict(row) for row in connection.execute(sql(name), parameters or {}).mappings()]


def json_default(value):
    if isinstance(value, Decimal):
        return float(value)
    if isinstance(value, (datetime, Path)):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=json_default) + "\n")


def import_regions(engine, records):
    """Commit the two imports together; a failed import leaves neither half done."""
    with engine.begin() as connection:
        present = connection.execute(text(
            "SELECT to_regclass('bioseq_external.sequences') IS NOT NULL"
        )).scalar_one()
        if not present:
            raise RuntimeError(
                "This BioseqDB build has no NCBI sequence tables. See the dated "
                "availability notice in README.md; a runtime SET cannot add them."
            )
        connection.execute(sql("02_schema.sql"))
        for record in records:
            connection.execute(text(
                "DELETE FROM butterflies.regions WHERE key = :key"
            ), {"key": record["key"]})
            connection.execute(sql("03_import.sql"), record)
            count = connection.execute(text(
                "SELECT count(*) FROM butterflies.regions WHERE key = :key"
            ), {"key": record["key"]}).scalar_one()
            if count != 1:
                raise RuntimeError(f"Expected one NCBI interval for {record['accession']}; got {count}")
            # Keep sequential requests comfortably below the public E-utilities
            # limit even if an upstream response is very fast.
            time.sleep(0.4)
        # Validate before leaving engine.begin(): bad ranges or changed bases
        # must roll back both replacements, preserving any previous valid data.
        validate_sequences(connection, records)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def validate_sequences(connection, records):
    stored = [dict(row) for row in connection.execute(text(
        "SELECT key, accession, start0, end0, seq::text AS sequence "
        "FROM butterflies.regions ORDER BY key"
    )).mappings()]
    expected = {record["key"]: record for record in records}
    selected = [row for row in stored if row["key"] in expected]
    if len(selected) != len(expected):
        raise AssertionError("Missing or duplicated tutorial intervals")
    sequences = {}
    for row in selected:
        record = expected[row["key"]]
        sequence = row["sequence"].upper()
        for field in ("accession", "start0", "end0"):
            if row[field] != record[field]:
                raise AssertionError(f"Stored {field} does not match sources.json")
        if len(sequence) != record["end0"] - record["start0"]:
            raise AssertionError("NCBI interval length does not match its coordinates")
        if set(sequence) - set("ACGTN"):
            raise AssertionError("Unexpected alphabet in the DNA5 example")
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if digest != record["sequence_sha256_uppercase"]:
            raise AssertionError(f"Sequence changed for pinned accession {record['accession']}")
        sequences[row["accession"]] = sequence
    return sequences


def validate_alignments(matches, sequences):
    """Check every coordinate and edit count against the actual fetched bases."""
    if not matches:
        raise AssertionError("Alignment returned no rows; refusing to draw empty charts")
    for row in matches:
        ts, te = row["target_start"], row["target_end"]
        qs, qe = row["query_start"], row["query_end"]
        target_full = sequences[row["target_accession"]]
        query_full = sequences[row["query_accession"]]
        if not (0 <= ts < te <= len(target_full) and 0 <= qs < qe <= len(query_full)):
            raise AssertionError("Alignment coordinates exceed the imported intervals")
        target = target_full[ts:te]
        query = query_full[qs:qe]
        if row["is_reverse"]:
            query = reverse_complement(query)
        operations = re.findall(r"(\d+)([MID=X])", row["cigar"])
        if not operations or "".join(n + op for n, op in operations) != row["cigar"]:
            raise AssertionError(f"Unexpected CIGAR syntax: {row['cigar']}")
        tpos = qpos = edits = columns = 0
        for raw_count, operation in operations:
            count = int(raw_count)
            columns += count
            if operation in "M=X":
                left, right = target[tpos:tpos + count], query[qpos:qpos + count]
                if len(left) != count or len(right) != count:
                    raise AssertionError("CIGAR consumes beyond an alignment span")
                mismatches = sum(a != b for a, b in zip(left, right))
                if (operation == "=" and mismatches) or (operation == "X" and mismatches != count):
                    raise AssertionError("CIGAR equality/mismatch operation disagrees with bases")
                edits += mismatches
                tpos += count
                qpos += count
            elif operation == "I":
                qpos += count
                edits += count
            else:
                tpos += count
                edits += count
        if (tpos, qpos, edits, columns) != (
            len(target), len(query), row["edits"], row["block_length"]
        ):
            raise AssertionError("CIGAR lengths or exact edit count disagree with BioseqDB")
        row["identity_percent"] = 100.0 * (columns - edits) / columns
    return {"checked_alignments": len(matches), "failures": 0}


def covered_bases(matches, prefix):
    intervals = sorted((row[f"{prefix}_start"], row[f"{prefix}_end"]) for row in matches)
    end = total = 0
    for start, stop in intervals:
        total += max(0, stop - max(start, end))
        end = max(end, stop)
    return total


def analyze(engine, records):
    with engine.connect() as connection:
        connection.execute(text("SET statement_timeout = '180s'"))
        sequences = validate_sequences(connection, records)
        preview = rows(connection, "04_inspect.sql")
        composition = rows(connection, "05_composition.sql")
        for row in composition:
            if sum(row[base] for base in "acgtn") != row["length"]:
                raise AssertionError("Base counts do not sum to the sequence length")
            record = next(record for record in records if record["accession"] == row["accession"])
            row["source_length"] = record["source_sequence_length_bp"]
            row["source_role"] = record["source_sequence_role"]
        started = time.perf_counter()
        matches = rows(connection, "06_align.sql")
        elapsed = time.perf_counter() - started
        gc_windows = rows(connection, "07_gc_windows.sql", {"window_size": 10000})
    validation = validate_alignments(matches, sequences)
    from .annotations import coding_sequence, compare_coding_sequences
    coding = {record["key"]: coding_sequence(sequences[record["accession"]], record)
              for record in records}
    coding_comparison = compare_coding_sequences(coding["erato"], coding["melpomene"])
    plotted = [row for row in matches if row["block_length"] >= 300
               and row["mapq"] >= 20 and row["identity_percent"] >= 70]
    if not plotted:
        raise AssertionError("No matches pass the stated plotting filter")
    summary = {
        "alignment_rows": len(matches),
        "reverse_rows": sum(row["is_reverse"] for row in matches),
        "plotted_rows": len(plotted),
        "plotted_reverse_rows": sum(row["is_reverse"] for row in plotted),
        "plot_filter": {"minimum_block_length": 300, "minimum_mapq": 20, "minimum_identity_percent": 70},
        "alignment_seconds_this_run": elapsed,
        "plotted_target_covered_bases": covered_bases(plotted, "target"),
        "plotted_query_covered_bases": covered_bases(plotted, "query"),
        "cigar_validation": validation,
        "coding_sequence_comparison": coding_comparison,
        "gc_window_size": 10000,
        "gc_window_count": len(gc_windows),
    }
    return {"preview": preview, "sequences": composition, "matches": matches,
            "plotted_matches": plotted, "gc_windows": gc_windows, "summary": summary}


def run(engine, output, reuse=False, render=True):
    manifest = json.loads((ROOT / "sources.json").read_text())
    records = manifest["records"]
    if not reuse:
        import_regions(engine, records)
    result = analyze(engine, records)
    from .annotations import load_annotations, load_optix_exons
    result["genes"] = load_annotations(records, output / "downloads")
    result["optix_exons"] = load_optix_exons(records, output / "downloads")
    expected_coding = manifest["coding_sequence_comparison"]
    for key, value in result["summary"]["coding_sequence_comparison"].items():
        if expected_coding[key] != value:
            raise AssertionError(f"Coding sequence validation changed: {key}")
    result["summary"]["analyzed_at"] = datetime.now(timezone.utc).isoformat()
    result["summary"]["source_retrieval_date"] = manifest["retrieval_date"]
    result["summary"]["versions"] = {
        name: importlib.metadata.version(name)
        for name in ("sqlalchemy", "psycopg", "plotly", "kaleido")
    }
    save_json(output / "analysis.json", result)
    save_json(output / "summary.json", result["summary"])
    if render:
        from .plots import write_figures
        write_figures(result["sequences"], result["plotted_matches"],
                      result["gc_windows"], result["genes"], output / "figures",
                      optix_exons=result["optix_exons"])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse", action="store_true", help="Recheck stored sequence hashes and rerun analysis without importing again")
    parser.add_argument("--no-figures", action="store_true", help="Run all SQL and validation without browser-based static export")
    parser.add_argument("--output", type=Path, default=ROOT / "output")
    args = parser.parse_args()
    url = os.environ.get("DATABASE_URL")
    if not url:
        parser.error("Export DATABASE_URL first (see .env.example)")
    # The current external provider needs custom plans to push bound accession
    # and interval values into NCBI requests. Keep psycopg auto-prepare disabled.
    engine = create_engine(url, pool_pre_ping=True,
                           connect_args={"prepare_threshold": None})
    try:
        result = run(engine, args.output, reuse=args.reuse, render=not args.no_figures)
        print(json.dumps(result["summary"], indent=2, default=json_default))
    finally:
        engine.dispose()


if __name__ == "__main__":
    main()
