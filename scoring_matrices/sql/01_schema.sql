\set ON_ERROR_STOP on
CREATE EXTENSION IF NOT EXISTS bioseq_cache;
CREATE EXTENSION IF NOT EXISTS bioseq;
CREATE SCHEMA IF NOT EXISTS scoring;
CREATE TABLE IF NOT EXISTS scoring.genomes (
    accession text NOT NULL,
    role text NOT NULL,
    seq bioseq.seq('dna5', 64) NOT NULL
) USING bioseq;

CREATE TABLE IF NOT EXISTS scoring.samples (
    accession text PRIMARY KEY,isolate text UNIQUE NOT NULL,
    organism text NOT NULL,bases integer NOT NULL,title text NOT NULL
);
TRUNCATE scoring.samples;
\copy scoring.samples FROM 'samples.csv' CSV HEADER
