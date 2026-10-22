\set ON_ERROR_STOP on
BEGIN;
-- The known projection contains an insertion, mismatch, deletion and refskip.
CREATE TEMP TABLE test_runs AS
WITH parsed AS (
 SELECT r.parts[1]::integer AS span,r.parts[2] AS op,r.i
 FROM regexp_matches('2=1I1=1X1D1N2M','([0-9]+)([MIDNSHP=X])','g')
 WITH ORDINALITY r(parts,i)
)
SELECT *,coalesce(sum(span) FILTER(WHERE op IN ('M','=','X','D','N'))
 OVER(ORDER BY i ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING),0)::integer AS ref0,
 coalesce(sum(span) FILTER(WHERE op IN ('M','=','X','I'))
 OVER(ORDER BY i ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING),0)::integer AS query0
FROM parsed;
DO $test$
DECLARE projected text;
BEGIN
 SELECT string_agg(CASE WHEN op IN ('M','=','X') THEN substr('ACAGCGT',query0+1,span)
    WHEN op IN ('D','N') THEN repeat('-',span) ELSE '' END,'' ORDER BY i)
 INTO projected FROM test_runs;
 IF projected<>'ACGC--GT' THEN RAISE EXCEPTION 'Known CIGAR projection failed: %',projected; END IF;
 IF (SELECT array_agg(ref0 ORDER BY i) FROM test_runs)<>ARRAY[0,2,2,3,4,5,6]
 OR (SELECT array_agg(query0 ORDER BY i) FROM test_runs)<>ARRAY[0,2,3,4,5,5,5] THEN
 RAISE EXCEPTION 'Known CIGAR coordinate accounting failed'; END IF;
END
$test$;
CREATE TEMP TABLE test_pairs AS
WITH letters AS (SELECT x AS base,count(*)::bigint AS n
 FROM unnest(ARRAY['A','A','C','N','-']) x WHERE x IN ('A','C','G','T') GROUP BY x)
SELECT a.base AS a,b.base AS b,
 CASE WHEN a.base=b.base THEN a.n*(a.n-1) ELSE a.n*b.n END AS n
FROM letters a CROSS JOIN letters b;
DO $test$
BEGIN
 IF (SELECT sum(n) FROM test_pairs)<>6
 OR (SELECT n FROM test_pairs WHERE a='A' AND b='A')<>2
 OR (SELECT n FROM test_pairs WHERE a='A' AND b='C')<>2
 OR EXISTS(SELECT 1 FROM test_pairs WHERE a='N' OR b='N' OR a='-' OR b='-') THEN
 RAISE EXCEPTION 'Known pair-count or unknown-base exclusion check failed'; END IF;
 IF NOT EXISTS(SELECT 1 FROM scoring.mitochondrial_code WHERE codon='AAT' AND aa='N')
 OR NOT EXISTS(SELECT 1 FROM scoring.mitochondrial_code WHERE codon='TGA' AND aa='W')
 OR NOT EXISTS(SELECT 1 FROM scoring.mitochondrial_code WHERE codon='ATA' AND aa='M') THEN
 RAISE EXCEPTION 'Mitochondrial translation/asparagine check failed'; END IF;
END
$test$;
ROLLBACK;
SELECT 'Toy CIGAR, pair-count and mitochondrial-code invariants passed' AS validation;
