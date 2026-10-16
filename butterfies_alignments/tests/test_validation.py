"""Regression checks for provenance and atomic imports.

Run ``poetry run python -m unittest discover -s tests -v``. Set DATABASE_URL to
also exercise rollback in a disposable PostgreSQL schema. The database test
uses tiny synthetic sequences as its remote-source fixture, never NCBI data.
"""
from copy import deepcopy
import gzip
import hashlib
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from uuid import uuid4

from sqlalchemy import create_engine, event, text

from butterfly_tutorial.annotations import load_annotations, load_optix_exons
from butterfly_tutorial.workflow import import_regions


class AnnotationValidation(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.cache = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)
        self.record = {
            "assembly": "GCF_TEST.1", "accession": "NC_TEST.1",
            "start0": 0, "end0": 100,
            "optix_gene": {"id": "gene-Optix", "start0": 9, "end0": 99, "strand": "+"},
            "optix_exons": [{"start0": 9, "end0": 39, "strand": "+"},
                            {"start0": 59, "end0": 99, "strand": "+"}],
            "optix_cds": [{"start0": 19, "end0": 28, "strand": "+", "protein_id": "XP_TEST.1"}],
        }
        self.gff = (
            "##gff-version 3\n"
            "NC_TEST.1\tTest\tgene\t10\t99\t.\t+\t.\tID=gene-Optix;Name=Optix\n"
            "NC_TEST.1\tTest\tmRNA\t10\t99\t.\t+\t.\tID=rna-Optix;Parent=gene-Optix\n"
            "NC_TEST.1\tTest\texon\t10\t39\t.\t+\t.\tID=exon-1;Parent=rna-Optix\n"
            "NC_TEST.1\tTest\texon\t60\t99\t.\t+\t.\tID=exon-2;Parent=rna-Optix\n"
            "NC_TEST.1\tTest\tCDS\t20\t28\t.\t+\t0\tID=cds-1;Parent=rna-Optix;protein_id=XP_TEST.1\n"
        )

    def write_gff(self, source=None):
        data = gzip.compress((self.gff if source is None else source).encode(), mtime=0)
        (self.cache / "GCF_TEST.1_genomic.gff.gz").write_bytes(data)
        # Updating this checksum deliberately tests relationships, not just the
        # easier compressed-file corruption check.
        self.record["gff_sha256"] = hashlib.sha256(data).hexdigest()

    def test_valid_gene_transcript_cds_and_exons(self):
        self.write_gff()
        self.assertEqual(len(load_annotations([self.record], self.cache)), 1)
        exons = load_optix_exons([self.record], self.cache)
        self.assertEqual([(e["start"], e["end"]) for e in exons], [(9, 39), (59, 99)])

    def test_cds_cannot_belong_to_another_transcript(self):
        self.write_gff(self.gff.replace("ID=cds-1;Parent=rna-Optix", "ID=cds-1;Parent=rna-Other"))
        with self.assertRaisesRegex(ValueError, "CDS features disagree"):
            load_annotations([self.record], self.cache)

    def test_exon_cannot_belong_to_another_transcript(self):
        self.write_gff(self.gff.replace("ID=exon-2;Parent=rna-Optix", "ID=exon-2;Parent=rna-Other"))
        with self.assertRaisesRegex(ValueError, "exon features disagree"):
            load_optix_exons([self.record], self.cache)

    def test_coordinate_mismatch_in_each_manifest_feature(self):
        self.write_gff()
        for feature in ("optix_exons", "optix_cds"):
            with self.subTest(feature=feature):
                record = deepcopy(self.record)
                record[feature][0]["start0"] += 1
                with self.assertRaisesRegex(ValueError, "features disagree"):
                    load_annotations([record], self.cache)

    def test_protein_accession_and_strand_are_checked(self):
        self.write_gff()
        for field, wrong in (("protein_id", "XP_OTHER.1"), ("strand", "-")):
            with self.subTest(field=field):
                record = deepcopy(self.record)
                record["optix_cds"][0][field] = wrong
                with self.assertRaisesRegex(ValueError, "CDS features disagree"):
                    load_annotations([record], self.cache)

    def test_unsupported_translation_models_are_rejected(self):
        for changed in (
            self.gff.replace("\t+\t0\tID=cds", "\t+\t1\tID=cds"),
            self.gff.replace("protein_id=XP_TEST.1", "protein_id=XP_TEST.1;partial=true"),
            self.gff.replace("protein_id=XP_TEST.1", "protein_id=XP_TEST.1;transl_except=test"),
            self.gff + "NC_TEST.1\tTest\tmRNA\t10\t99\t.\t+\t.\tID=rna-Second;Parent=gene-Optix\n",
        ):
            with self.subTest(gff=changed):
                self.write_gff(changed)
                with self.assertRaises(ValueError):
                    load_annotations([self.record], self.cache)


@unittest.skipUnless(os.environ.get("DATABASE_URL"), "Set DATABASE_URL for the PostgreSQL rollback regression")
class ImportRollback(unittest.TestCase):
    def setUp(self):
        self.schema = "butterfly_test_" + uuid4().hex
        self.engine = create_engine(os.environ["DATABASE_URL"], connect_args={"prepare_threshold": None})
        self.addCleanup(self.cleanup_database)
        with self.engine.begin() as connection:
            connection.execute(text(f"CREATE SCHEMA {self.schema}"))
            connection.execute(text(f"""CREATE TABLE {self.schema}.regions (
                key text, accession text, organism text, assembly text,
                start0 bigint, end0 bigint, seq bioseq.seq('dna5', 64))"""))
            connection.execute(text(f"""CREATE TABLE {self.schema}.remote_sequences (
                provider text, accession text, start0 bigint, end0 bigint,
                sequence bioseq.seq('dna5', 64))"""))
            for i, sequence in enumerate(("ACGT", "GCAA")):
                connection.execute(text(f"""INSERT INTO {self.schema}.regions VALUES
                    (:key, :accession, 'Synthetic organism', 'TEST.1', 0, 4,
                     'TTTT'::bioseq.seq('dna5', 64))"""), {"key": f"test{i}", "accession": f"TEST{i}.1"})
                connection.execute(text(f"""INSERT INTO {self.schema}.remote_sequences VALUES
                    ('ncbi', :accession, 0, 4, CAST(:sequence AS bioseq.seq('dna5', 64)))"""),
                    {"accession": f"TEST{i}.1", "sequence": sequence})
        self.records = [{"key": f"test{i}", "accession": f"TEST{i}.1", "organism": "Synthetic organism",
                         "assembly": "TEST.1", "start0": 0, "end0": 4,
                         "sequence_sha256_uppercase": hashlib.sha256(sequence.encode()).hexdigest()}
                        for i, sequence in enumerate(("ACGT", "GCAA"))]
        # Only this test engine rewrites the two fixed tutorial relations. No
        # operation can reach the real butterflies schema or a remote service.
        self.rewrite_listener = self.rewrite_relations
        event.listen(self.engine, "before_cursor_execute", self.rewrite_listener, retval=True)

    def rewrite_relations(self, connection, cursor, statement, parameters, context, executemany):
        statement = statement.replace("bioseq_external.sequences", f"{self.schema}.remote_sequences")
        statement = statement.replace("butterflies.", f"{self.schema}.")
        statement = statement.replace("CREATE SCHEMA IF NOT EXISTS butterflies", f"CREATE SCHEMA IF NOT EXISTS {self.schema}")
        return statement, parameters

    def cleanup_database(self):
        if hasattr(self, "rewrite_listener"):
            event.remove(self.engine, "before_cursor_execute", self.rewrite_listener)
        with self.engine.begin() as connection:
            connection.execute(text(f"DROP SCHEMA IF EXISTS {self.schema} CASCADE"))
        self.engine.dispose()

    def sequences(self):
        with self.engine.connect() as connection:
            return list(connection.execute(text(f"SELECT key, seq::text FROM {self.schema}.regions ORDER BY key")))

    def test_hash_mismatch_rolls_back_both_replacements(self):
        bad_records = deepcopy(self.records)
        bad_records[1]["sequence_sha256_uppercase"] = "0" * 64
        with patch("butterfly_tutorial.workflow.time.sleep"), self.assertRaisesRegex(AssertionError, "Sequence changed"):
            import_regions(self.engine, bad_records)
        self.assertEqual(self.sequences(), [("test0", "TTTT"), ("test1", "TTTT")])
        # A valid run must really replace both rows; otherwise an implementation
        # that never imports anything could accidentally satisfy the rollback test.
        with patch("butterfly_tutorial.workflow.time.sleep"):
            import_regions(self.engine, self.records)
        self.assertEqual(self.sequences(), [("test0", "ACGT"), ("test1", "GCAA")])


if __name__ == "__main__":
    unittest.main()
