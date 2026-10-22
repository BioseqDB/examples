\set ON_ERROR_STOP on
-- Recalculate edit counts directly from aligned letters and CIGAR insertions.
-- BioseqDB counts identical N/N letters as equal for its edit field.
-- All N letters are independently excluded from biological pair counts.
DROP TABLE IF EXISTS scoring.cigar_validation;
CREATE TABLE scoring.cigar_validation AS
WITH letter_edits AS (
    SELECT p.accession,count(*) FILTER(
        WHERE b.base='-' OR b.base<>r.base) AS changes
    FROM scoring.projected_sequences p
    CROSS JOIN LATERAL regexp_split_to_table(p.aligned,'') WITH ORDINALITY b(base,pos)
    JOIN scoring.reference_bases r ON r.position1=b.pos
    JOIN scoring.best_alignment a ON a.accession=p.accession
    WHERE b.pos>a.target_start AND b.pos<=a.target_end
    GROUP BY p.accession
), insertions AS (
    SELECT accession,coalesce(sum(span) FILTER(WHERE op='I'),0) AS bases
    FROM scoring.cigar_runs GROUP BY accession
)
SELECT a.accession,a.edits,l.changes+i.bases AS replayed_edits,
       (SELECT sum(span) FROM scoring.cigar_runs r WHERE r.accession=a.accession
          AND op IN ('M','=','X','D','N'))=a.target_end-a.target_start AS reference_span_ok,
       (SELECT sum(span) FROM scoring.cigar_runs r WHERE r.accession=a.accession
          AND op IN ('M','=','X','I'))=a.query_end-a.query_start AS query_span_ok
FROM scoring.best_alignment a JOIN letter_edits l USING(accession)
JOIN insertions i USING(accession);

-- Independently translate the reference annotation; compare all 13 proteins.
DROP TABLE IF EXISTS scoring.reference_translation_check;
CREATE TABLE scoring.reference_translation_check AS
WITH coding AS (
    SELECT c.gene,c.reference_protein,
       CASE WHEN c.strand='+' THEN substr(g.seq::text,c.start0+1,c.end0-c.start0)
       ELSE reverse(translate(substr(g.seq::text,c.start0+1,c.end0-c.start0),'ACGT','TGCA')) END AS bases
    FROM scoring.cds c CROSS JOIN scoring.genomes g WHERE g.role='reference'
), codons AS (
    SELECT c.*,i,substr(bases,3*i-2,3) AS codon
    FROM coding c CROSS JOIN LATERAL generate_series(1,length(reference_protein)) x(i)
), translated AS (
    SELECT gene,reference_protein,string_agg(
      CASE WHEN i=1 AND codon IN ('ATT','ATC','ATA','ATG','GTG') THEN 'M'
           ELSE coalesce(t.aa,'X') END,'' ORDER BY i) AS computed_protein
    FROM codons LEFT JOIN scoring.mitochondrial_code t USING(codon)
    GROUP BY gene,reference_protein
)
SELECT gene,length(reference_protein) AS amino_acids,
       computed_protein=reference_protein AS exact_match FROM translated;

DO $validate$
BEGIN
    IF (SELECT count(*) FROM scoring.genomes WHERE role='sample')<>1000 THEN
        RAISE EXCEPTION 'Expected exactly 1,000 sample genomes';
    END IF;
    IF (SELECT count(DISTINCT isolate) FROM scoring.samples)<>1000 THEN
        RAISE EXCEPTION 'Expected 1,000 distinct submitted sample isolates';
    END IF;
    IF EXISTS (SELECT accession FROM scoring.genomes WHERE role='sample'
               EXCEPT SELECT accession FROM scoring.samples)
       OR EXISTS (SELECT accession FROM scoring.samples
                  EXCEPT SELECT accession FROM scoring.genomes WHERE role='sample') THEN
        RAISE EXCEPTION 'Downloaded cohort differs from the accession manifest';
    END IF;
    IF EXISTS (SELECT 1 FROM scoring.genomes WHERE seq::text !~ '^[ACGTN]+$') THEN
        RAISE EXCEPTION 'Unexpected nucleotide alphabet';
    END IF;
    IF EXISTS (SELECT 1 FROM scoring.genomes g JOIN scoring.samples s USING(accession)
               WHERE bioseq.len(g.seq)<>s.bases OR s.organism<>'Homo sapiens'
                  OR s.title NOT LIKE '%complete genome%') THEN
        RAISE EXCEPTION 'Downloaded sequence length or metadata mismatch';
    END IF;
    IF EXISTS (SELECT accession FROM scoring.genomes GROUP BY accession HAVING count(*)<>1) THEN
        RAISE EXCEPTION 'Duplicate accession';
    END IF;
    IF (SELECT count(*) FROM scoring.best_alignment)<>1000 THEN
        RAISE EXCEPTION 'At least one sample lacks a forward primary alignment';
    END IF;
    IF EXISTS (SELECT 1 FROM scoring.best_alignment a JOIN scoring.genomes g USING(accession)
               WHERE target_end-target_start < 0.999*16569 OR query_start<>0
                  OR query_end<>bioseq.len(g.seq)) THEN
        RAISE EXCEPTION 'Expected full query coverage and at least 99.9 percent reference coverage';
    END IF;
    IF EXISTS (SELECT 1 FROM scoring.cigar_runs WHERE op NOT IN ('M','=','X','I','D')) THEN
        RAISE EXCEPTION 'Unexpected CIGAR operation';
    END IF;
    IF EXISTS (SELECT 1 FROM scoring.projected_sequences WHERE length(aligned)<>16569) THEN
        RAISE EXCEPTION 'Incorrect reference projection length';
    END IF;
    IF EXISTS (SELECT 1 FROM scoring.cigar_validation WHERE edits<>replayed_edits
               OR NOT reference_span_ok OR NOT query_span_ok) THEN
        RAISE EXCEPTION 'CIGAR replay mismatch';
    END IF;
    IF EXISTS (SELECT 1 FROM scoring.reference_translation_check WHERE NOT exact_match)
       OR (SELECT count(*) FROM scoring.reference_translation_check)<>13 THEN
        RAISE EXCEPTION 'Reference translation disagrees with NCBI annotation';
    END IF;
    IF (SELECT sum(n) FROM scoring.site_counts)<>16569000 THEN
        RAISE EXCEPTION 'Aligned column count mismatch';
    END IF;
    IF (SELECT sum(pairs) FROM scoring.nucleotide_pairs) <>
       (SELECT sum(n*(n-1)) FROM (
          SELECT sum(s.n) AS n FROM scoring.site_counts s
          JOIN scoring.reference_bases r USING(position1)
          WHERE s.base IN ('A','C','G','T') AND r.base IN ('A','C','G','T')
          GROUP BY s.position1) counts)
       OR (SELECT sum(pairs) FROM scoring.protein_pairs) <>
       (SELECT sum(n*(n-1)) FROM (
          SELECT sum(n) AS n FROM scoring.protein_site_counts GROUP BY gene,codon1) counts) THEN
        RAISE EXCEPTION 'Pair totals disagree with per-column callable denominators';
    END IF;
    IF abs((SELECT sum(q) FROM scoring.nucleotide_matrix)-1)>1e-12
       OR abs((SELECT sum(q) FROM scoring.protein_matrix)-1)>1e-12 THEN
        RAISE EXCEPTION 'Matrix frequencies do not sum to 1';
    END IF;
    IF EXISTS (SELECT 1 FROM scoring.nucleotide_pairs a JOIN scoring.nucleotide_pairs b
               ON a.base_from=b.base_to AND a.base_to=b.base_from WHERE a.pairs<>b.pairs)
       OR EXISTS (SELECT 1 FROM scoring.protein_pairs a JOIN scoring.protein_pairs b
               ON a.aa_from=b.aa_to AND a.aa_to=b.aa_from WHERE a.pairs<>b.pairs) THEN
        RAISE EXCEPTION 'Pair matrix is not symmetric';
    END IF;
    IF NOT EXISTS (SELECT 1 FROM scoring.nucleotide_pairs WHERE base_from<>base_to AND pairs>0)
       OR NOT EXISTS (SELECT 1 FROM scoring.protein_pairs WHERE aa_from<>aa_to AND pairs>0) THEN
        RAISE EXCEPTION 'No substitution observations';
    END IF;
END
$validate$;
SELECT 'All cohort, alignment, translation and matrix checks passed' AS validation;
