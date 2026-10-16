INSERT INTO butterflies.regions
    (key, accession, organism, assembly, start0, end0, seq)
SELECT :key, :accession, :organism, :assembly, :start0, :end0,
       sequence::bioseq.seq('dna5', 64)
FROM bioseq_external.sequences
WHERE provider = 'ncbi'
  AND accession = :accession
  AND start0 = :start0
  AND end0 = :end0;
