\set ON_ERROR_STOP on
CREATE TABLE IF NOT EXISTS scoring.reference_joint (
    matrix text,aa_from text,aa_to text,q numeric,PRIMARY KEY(matrix,aa_from,aa_to)
);
CREATE TABLE IF NOT EXISTS scoring.reference_scores (
    matrix text,aa_from text,aa_to text,score integer,units_per_bit integer,
    PRIMARY KEY(matrix,aa_from,aa_to)
);
TRUNCATE scoring.reference_joint,scoring.reference_scores;
CREATE TEMP TABLE incoming_joint(aa_from text,aa_to text,q numeric);
CREATE TEMP TABLE incoming_scores(aa_from text,aa_to text,score integer,units_per_bit integer);
\copy incoming_joint FROM 'reference_matrices/blosum62_joint_probabilities.csv' CSV HEADER
INSERT INTO scoring.reference_joint SELECT 'BLOSUM62',* FROM incoming_joint;
TRUNCATE incoming_joint;
\copy incoming_joint FROM 'reference_matrices/pam250_joint_probabilities.csv' CSV HEADER
INSERT INTO scoring.reference_joint SELECT 'PAM250',* FROM incoming_joint;
\copy incoming_scores FROM 'reference_matrices/blosum62_scores.csv' CSV HEADER
INSERT INTO scoring.reference_scores SELECT 'BLOSUM62',* FROM incoming_scores;
TRUNCATE incoming_scores;
\copy incoming_scores FROM 'reference_matrices/pam250_scores.csv' CSV HEADER
INSERT INTO scoring.reference_scores SELECT 'PAM250',* FROM incoming_scores;

CREATE OR REPLACE VIEW scoring.reference_protein_matrix AS
WITH background AS (
    SELECT matrix,aa_from AS aa,sum(q) AS p
    FROM scoring.reference_joint GROUP BY matrix,aa_from
)
SELECT j.matrix,j.aa_from,j.aa_to,j.q,a.p*b.p AS expected,
       s.score,s.units_per_bit,2*log(2::numeric,j.q/(a.p*b.p)) AS score_half_bits,
       CASE WHEN j.aa_from<>j.aa_to THEN
           j.q/nullif(sum(j.q) FILTER(WHERE j.aa_from<>j.aa_to)
             OVER(PARTITION BY j.matrix,j.aa_from),0)
       END AS fraction_of_differences
FROM scoring.reference_joint j JOIN scoring.reference_scores s USING(matrix,aa_from,aa_to)
JOIN background a ON a.matrix=j.matrix AND a.aa=j.aa_from
JOIN background b ON b.matrix=j.matrix AND b.aa=j.aa_to;

CREATE OR REPLACE VIEW scoring.matrix_entropy AS
SELECT 'Human mitochondrial DNA'::text AS matrix,
       sum(q*score_half_bits/2) AS bits_per_aligned_pair FROM scoring.nucleotide_matrix
UNION ALL
SELECT 'Human mitochondrial proteins',sum(q*score_half_bits/2) FROM scoring.protein_matrix
UNION ALL
SELECT matrix,sum(q*score_half_bits/2) FROM scoring.reference_protein_matrix GROUP BY matrix;
