\set ON_ERROR_STOP on
SET statement_timeout = '20min';
DROP TABLE IF EXISTS scoring.alignments CASCADE;
CREATE TABLE scoring.alignments AS
SELECT q.accession,
       bioseq.target_start(t.seq)::integer AS target_start,
       bioseq.target_end(t.seq)::integer AS target_end,
       bioseq.query_start(t.seq)::integer AS query_start,
       bioseq.query_end(t.seq)::integer AS query_end,
       bioseq.is_reverse(t.seq) AS is_reverse,
       bioseq.cigar(t.seq)::text AS cigar,
       bioseq.edit(t.seq) AS edits,
       bioseq.score(t.seq) AS score,
       bioseq.is_primary(t.seq) AS is_primary
FROM scoring.genomes t
JOIN scoring.genomes q ON mmsearch(t.seq, q.seq, aligner_options(
    kmer_size => 15, window_size => 10,
    max_gap => 200, min_chain_score => 40))
WHERE t.role = 'reference' AND q.role = 'sample';

-- Choose one forward primary alignment per genome. Any unexpected split,
-- reversed or uncovered sequence is visible in the quality report.
CREATE OR REPLACE VIEW scoring.best_alignment AS
SELECT DISTINCT ON (accession) *
FROM scoring.alignments
WHERE is_primary AND NOT is_reverse
ORDER BY accession, score DESC, target_start, query_start;
SELECT count(*) AS aligned_samples,
       min(target_end - target_start) AS shortest_reference_span,
       max(target_end - target_start) AS longest_reference_span,
       min(edits) AS fewest_edits, max(edits) AS most_edits
FROM scoring.best_alignment;
