\set ON_ERROR_STOP on
DROP TABLE IF EXISTS scoring.nucleotide_pairs CASCADE;
CREATE TABLE scoring.nucleotide_pairs AS
SELECT a.base AS base_from,b.base AS base_to,
       sum(CASE WHEN a.base=b.base THEN a.n*(a.n-1)
                ELSE a.n*b.n END)::bigint AS pairs
FROM scoring.site_counts a JOIN scoring.site_counts b USING(position1)
JOIN scoring.reference_bases r USING(position1)
WHERE r.base IN ('A','C','G','T') AND a.base IN ('A','C','G','T') AND b.base IN ('A','C','G','T')
GROUP BY a.base,b.base;

-- Ordered pairs give a symmetric matrix and use e_ij = p_i * p_j.
-- Two orientations represent each unordered pair of different samples.
CREATE OR REPLACE VIEW scoring.nucleotide_matrix AS
WITH frequencies AS (
    SELECT *, pairs::numeric/sum(pairs) OVER () AS q
    FROM scoring.nucleotide_pairs
), marginals AS (
    SELECT base_from AS base, sum(q) AS p FROM frequencies GROUP BY base_from
)
SELECT f.base_from, f.base_to, pairs,q,a.p*b.p AS expected,
       CASE WHEN q>0 THEN 2*log(2::numeric,q/(a.p*b.p)) END AS score_half_bits,
       CASE WHEN f.base_from<>f.base_to THEN
           pairs::numeric/nullif(sum(pairs) FILTER (WHERE f.base_from<>f.base_to)
               OVER (PARTITION BY f.base_from),0)
       END AS fraction_of_differences
FROM frequencies f JOIN marginals a ON a.base=f.base_from
JOIN marginals b ON b.base=f.base_to;

CREATE OR REPLACE VIEW scoring.reference_spectrum AS
WITH tallies AS (
    SELECT r.base AS reference_base,s.base AS sample_base,sum(s.n)::bigint AS calls
    FROM scoring.site_counts s JOIN scoring.reference_bases r USING(position1)
    WHERE r.base IN ('A','C','G','T') AND s.base IN ('A','C','G','T')
    GROUP BY r.base,s.base
)
SELECT *, calls::numeric/sum(calls) OVER(PARTITION BY reference_base) AS fraction_of_calls,
       CASE WHEN reference_base<>sample_base THEN
         calls::numeric/nullif(sum(calls) FILTER(WHERE reference_base<>sample_base)
             OVER(PARTITION BY reference_base),0)
       END AS fraction_of_differences
FROM tallies;

CREATE OR REPLACE VIEW scoring.position_summary AS
SELECT s.position1,r.base AS reference_base,
       sum(s.n) FILTER(WHERE s.base IN ('A','C','G','T')) AS called_samples,
       coalesce(sum(s.n) FILTER(WHERE s.base='N'),0) AS unknown_n,
       coalesce(sum(s.n) FILTER(WHERE s.base='-'),0) AS gaps,
       coalesce(sum(s.n) FILTER(WHERE r.base IN ('A','C','G','T') AND s.base IN ('A','C','G','T') AND s.base<>r.base),0) AS reference_differences,
       count(*) FILTER(WHERE s.base IN ('A','C','G','T')) AS observed_bases
FROM scoring.site_counts s JOIN scoring.reference_bases r USING(position1)
GROUP BY s.position1,r.base;

SELECT base_from,base_to,pairs,round(score_half_bits,2) AS score,
       round(100*fraction_of_differences,2) AS difference_percent
FROM scoring.nucleotide_matrix ORDER BY base_from,base_to;
