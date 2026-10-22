\set ON_ERROR_STOP on
DO $check$
BEGIN
    IF (SELECT count(*) FROM scoring.genomes WHERE role='sample')<>1000
       OR (SELECT count(DISTINCT accession) FROM scoring.genomes WHERE role='sample')<>1000
       OR (SELECT count(DISTINCT isolate) FROM scoring.samples)<>1000 THEN
        RAISE EXCEPTION 'Import is incomplete: expected 1,000 distinct sample genomes and isolates';
    END IF;
    IF (SELECT count(*) FROM scoring.genomes WHERE role='reference')<>1 THEN
        RAISE EXCEPTION 'Expected one reference genome';
    END IF;
END
$check$;
