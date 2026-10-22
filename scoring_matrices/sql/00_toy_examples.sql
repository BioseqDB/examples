\set ON_ERROR_STOP on

-- Additive log-odds scores
SELECT round(2 * log(2::numeric, 0.04 / 0.01)) AS half_bit_score;

SELECT sum(score) AS alignment_score
FROM (VALUES ('A', 'A', 4), ('L', 'I', 2), ('K', 'D', -1))
     AS aligned(first_residue, second_residue, score);

-- Unordered pairs in one aligned column
WITH residues(letter) AS (
    VALUES ('A'), ('A'), ('A'), ('A'), ('I'), ('L')
), counts AS (
    SELECT letter, count(*) AS n FROM residues GROUP BY letter
), pairs AS (
    SELECT a.letter AS first_residue, b.letter AS second_residue,
           CASE WHEN a.letter = b.letter THEN a.n * (a.n - 1) / 2
                ELSE a.n * b.n END AS observations
    FROM counts a JOIN counts b ON a.letter <= b.letter
)
SELECT * FROM pairs WHERE observations > 0
ORDER BY first_residue, second_residue;
