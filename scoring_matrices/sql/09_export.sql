\set ON_ERROR_STOP on
CREATE OR REPLACE VIEW scoring.genome_quality AS
WITH bases AS MATERIALIZED (SELECT accession,role,seq::text AS bases FROM scoring.genomes)
SELECT accession,role,length(bases) AS bp,
       length(bases)-length(replace(bases,'N','')) AS n_count,
       length(bases)-length(translate(bases,'ACGT','')) AS called_bases,
       length(bases)-length(translate(bases,'GC','')) AS gc_bases,
       100.0*(length(bases)-length(translate(bases,'GC','')))/
         nullif(length(bases)-length(translate(bases,'ACGT','')),0) AS gc_percent
FROM bases;

CREATE OR REPLACE VIEW scoring.cohort_quality AS
SELECT role,count(*) AS genomes,min(bp) AS shortest,max(bp) AS longest,
       sum(bp) AS total_bases,sum(n_count) AS n_count,sum(called_bases) AS called_bases,
       100.0*sum(gc_bases)/nullif(sum(called_bases),0) AS gc_percent
FROM scoring.genome_quality GROUP BY role;

\copy (SELECT * FROM scoring.cohort_quality ORDER BY role) TO 'results/quality.csv' CSV HEADER
\copy (SELECT * FROM scoring.genome_quality ORDER BY accession) TO 'results/genome_quality.csv' CSV HEADER
\copy (SELECT * FROM scoring.best_alignment ORDER BY accession) TO 'results/alignment_quality.csv' CSV HEADER
\copy (SELECT * FROM scoring.nucleotide_matrix ORDER BY base_from,base_to) TO 'results/nucleotide_matrix.csv' CSV HEADER
\copy (SELECT * FROM scoring.reference_spectrum ORDER BY reference_base,sample_base) TO 'results/reference_spectrum.csv' CSV HEADER
\copy (SELECT * FROM scoring.position_summary ORDER BY position1) TO 'results/position_summary.csv' CSV HEADER
\copy (SELECT * FROM scoring.protein_matrix ORDER BY aa_from,aa_to) TO 'results/protein_matrix.csv' CSV HEADER
\copy (SELECT * FROM scoring.reference_protein_matrix ORDER BY matrix,aa_from,aa_to) TO 'results/reference_protein_matrix.csv' CSV HEADER
\copy (SELECT * FROM scoring.matrix_entropy ORDER BY matrix) TO 'results/matrix_entropy.csv' CSV HEADER
\copy (SELECT c.gene,count(*) AS sample_cds,count(*) FILTER(WHERE has_indel) AS indel_cds,count(*) FILTER(WHERE incompletely_aligned) AS incomplete_cds,count(*) FILTER(WHERE NOT has_indel AND NOT incompletely_aligned AND NOT q.premature_stop AND NOT q.invalid_start AND NOT q.ambiguous_start) AS retained_cds,count(*) FILTER(WHERE q.premature_stop) AS stop_cds,count(*) FILTER(WHERE q.invalid_start) AS invalid_start_cds,count(*) FILTER(WHERE q.ambiguous_start) AS ambiguous_start_cds FROM scoring.cds_quality c LEFT JOIN scoring.protein_cds_quality q USING(accession,gene) GROUP BY c.gene ORDER BY c.gene) TO 'results/cds_quality.csv' CSV HEADER
\copy (SELECT count(*) AS evaluated_codons,count(*) FILTER(WHERE aa IS NULL) AS ambiguous_codons,count(*) FILTER(WHERE aa='*') AS stop_codons,count(*) FILTER(WHERE aa IS NOT NULL AND aa<>'*' AND NOT q.premature_stop AND NOT q.invalid_start AND NOT q.ambiguous_start) AS callable_codons,count(*) FILTER(WHERE aa='N' AND NOT q.premature_stop AND NOT q.invalid_start AND NOT q.ambiguous_start) AS asparagine_codons,(SELECT count(*) FROM (SELECT gene,codon1 FROM scoring.protein_site_counts GROUP BY gene,codon1 HAVING count(*)>1) s) AS polymorphic_codons,(SELECT count(*) FROM scoring.protein_cds_quality WHERE premature_stop) AS excluded_stop_cds,(SELECT count(*) FROM scoring.protein_cds_quality WHERE invalid_start) AS excluded_invalid_start_cds,(SELECT count(*) FROM scoring.protein_cds_quality WHERE ambiguous_start) AS excluded_ambiguous_start_cds,count(*) FILTER(WHERE NOT q.premature_stop AND NOT q.invalid_start AND NOT q.ambiguous_start) AS eligible_codons,count(*) FILTER(WHERE aa IS NULL AND NOT q.premature_stop AND NOT q.invalid_start AND NOT q.ambiguous_start) AS ambiguous_eligible_codons,count(*) FILTER(WHERE q.premature_stop OR q.invalid_start OR q.ambiguous_start) AS excluded_cds_codons FROM scoring.amino_acids a JOIN scoring.protein_cds_quality q USING(accession,gene)) TO 'results/protein_quality.csv' CSV HEADER
\copy (SELECT * FROM scoring.sample_distances ORDER BY sample_from,sample_to) TO 'results/sample_distances.csv' CSV HEADER
\copy (SELECT base,sum(n) AS aligned_calls FROM scoring.site_counts GROUP BY base ORDER BY base) TO 'results/aligned_base_counts.csv' CSV HEADER
\copy (SELECT op,count(*) AS runs,sum(span) AS bases FROM scoring.cigar_runs GROUP BY op ORDER BY op) TO 'results/cigar_summary.csv' CSV HEADER
