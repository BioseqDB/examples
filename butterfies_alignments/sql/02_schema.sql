CREATE SCHEMA IF NOT EXISTS butterflies;

CREATE TABLE IF NOT EXISTS butterflies.regions (
    key text NOT NULL,
    accession text NOT NULL,
    organism text NOT NULL,
    assembly text NOT NULL,
    start0 bigint NOT NULL,
    end0 bigint NOT NULL,
    seq bioseq.seq('dna5', 64) NOT NULL
) USING bioseq;
