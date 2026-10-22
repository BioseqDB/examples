\set ON_ERROR_STOP on
SELECT role, count(*) AS genomes,
       min(bioseq.len(seq)) AS shortest,
       max(bioseq.len(seq)) AS longest,
       sum(bioseq.len(seq)) AS total_bases,
       sum(length(seq::text) - length(replace(seq::text, 'N', ''))) AS unknown_n
FROM scoring.genomes GROUP BY role ORDER BY role;

SELECT accession, bioseq.len(seq) AS bases,
       length(seq::text) - length(replace(seq::text, 'N', '')) AS unknown_n
FROM scoring.genomes WHERE role = 'sample'
ORDER BY accession LIMIT 5;
