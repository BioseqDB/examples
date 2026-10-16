SELECT accession, scientific_name, assembly_level, assembly_status
FROM bioseq_external.catalog
WHERE provider = 'ncbi'
  AND scientific_name = :organism
ORDER BY accession;
