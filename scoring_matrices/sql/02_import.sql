\set ON_ERROR_STOP on
-- Each accession is a versioned record from Li et al. (2014).
-- The first 1,000 complete-genome records form a reproducible tutorial subset.
-- Keep at least 0.4 seconds between requests; rerunning resumes missing rows.
DELETE FROM scoring.genomes g
WHERE g.role='sample' AND NOT EXISTS (
    SELECT 1 FROM scoring.samples s WHERE s.accession=g.accession
);

WITH requested AS (
    SELECT 'NC_012920.1' AS accession, 'reference' AS role
    UNION ALL
    SELECT accession, 'sample' FROM scoring.samples
)
SELECT format($fetch$
    INSERT INTO scoring.genomes (accession, role, seq)
    SELECT %L, %L, sequence::bioseq.seq('dna5', 64)
    FROM bioseq_external.sequences
    WHERE provider = 'ncbi' AND accession = %L;
    SELECT pg_sleep(0.4);
$fetch$, r.accession, r.role, r.accession)
FROM requested r
WHERE NOT EXISTS (
    SELECT 1 FROM scoring.genomes g WHERE g.accession = r.accession
)
ORDER BY r.role, r.accession
\gexec
