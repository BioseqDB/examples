\set ON_ERROR_STOP on
DROP TABLE IF EXISTS scoring.cigar_runs CASCADE;
CREATE TABLE scoring.cigar_runs AS
WITH operations AS (
    SELECT a.*, r.ordinality AS run,
           r.parts[1]::integer AS span, r.parts[2] AS op
    FROM scoring.best_alignment a
    CROSS JOIN LATERAL regexp_matches(a.cigar, '([0-9]+)([MIDNSHP=X])', 'g')
        WITH ORDINALITY AS r(parts, ordinality)
)
SELECT accession, run, span, op,
       target_start + coalesce(sum(span) FILTER (WHERE op IN ('M','=','X','D','N'))
         OVER (PARTITION BY accession ORDER BY run ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING),0)::integer AS ref_start0,
       query_start + coalesce(sum(span) FILTER (WHERE op IN ('M','=','X','I'))
         OVER (PARTITION BY accession ORDER BY run ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING),0)::integer AS query_start0
FROM operations;

DROP TABLE IF EXISTS scoring.projected_sequences CASCADE;
CREATE TABLE scoring.projected_sequences AS
WITH sequences AS MATERIALIZED (
    SELECT accession, seq::text AS bases FROM scoring.genomes WHERE role = 'sample'
), blocks AS (
    SELECT r.accession, string_agg(
        CASE WHEN op IN ('M','=','X') THEN substr(g.bases,r.query_start0+1,r.span)
             WHEN op IN ('D','N') THEN repeat('-',r.span)
             ELSE '' END, '' ORDER BY run) AS aligned
    FROM scoring.cigar_runs r JOIN sequences g USING (accession)
    GROUP BY r.accession
)
SELECT a.accession,
       repeat('-',a.target_start) || b.aligned ||
       repeat('-', bioseq.len(g.seq)::integer-a.target_end) AS aligned
FROM scoring.best_alignment a JOIN blocks b USING (accession)
CROSS JOIN scoring.genomes g WHERE g.role = 'reference';

DROP TABLE IF EXISTS scoring.reference_bases CASCADE;
CREATE TABLE scoring.reference_bases AS
SELECT p.pos::integer AS position1, p.base
FROM scoring.genomes g
CROSS JOIN LATERAL regexp_split_to_table(g.seq::text,'')
  WITH ORDINALITY AS p(base,pos)
WHERE g.role='reference';

DROP TABLE IF EXISTS scoring.site_counts CASCADE;
CREATE TABLE scoring.site_counts AS
SELECT p.pos::integer AS position1, p.base, count(*)::bigint AS n
FROM scoring.projected_sequences g
CROSS JOIN LATERAL regexp_split_to_table(g.aligned,'')
  WITH ORDINALITY AS p(base,pos)
GROUP BY p.pos,p.base;
CREATE UNIQUE INDEX ON scoring.site_counts(position1,base);

SELECT base, sum(n) AS aligned_calls
FROM scoring.site_counts GROUP BY base ORDER BY base;
