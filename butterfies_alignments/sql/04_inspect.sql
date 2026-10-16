SELECT key, accession, bioseq.len(seq) AS length,
       seq[0:20]::text AS first_20_bases
FROM butterflies.regions
ORDER BY key;
