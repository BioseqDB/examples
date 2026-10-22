\set ON_ERROR_STOP on
DROP TABLE IF EXISTS scoring.sample_distances;
CREATE TABLE scoring.sample_distances AS
WITH chosen AS MATERIALIZED (
    SELECT p.accession,r.pos::integer AS position1,r.base
    FROM scoring.projected_sequences p
    CROSS JOIN LATERAL regexp_split_to_table(p.aligned,'') WITH ORDINALITY r(base,pos)
    WHERE p.accession IN (SELECT 'KF'||(161060+142*i)||'.1' FROM generate_series(0,7) x(i))
), counts AS (
    SELECT a.accession AS sample_from,b.accession AS sample_to,
           count(*) AS compared_sites,count(*) FILTER(WHERE a.base<>b.base) AS differences
    FROM chosen a JOIN chosen b USING(position1)
    JOIN scoring.reference_bases r USING(position1)
    WHERE a.base IN ('A','C','G','T') AND b.base IN ('A','C','G','T')
      AND r.base IN ('A','C','G','T')
    GROUP BY a.accession,b.accession
)
SELECT *,differences::numeric/compared_sites AS p_distance FROM counts;
