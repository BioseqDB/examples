WITH counts AS (
    SELECT key, accession, organism, start0 AS interval_start,
           bioseq.len(seq) AS length,
           bioseq.symbol_count(seq, 'A') AS a,
           bioseq.symbol_count(seq, 'C') AS c,
           bioseq.symbol_count(seq, 'G') AS g,
           bioseq.symbol_count(seq, 'T') AS t,
           bioseq.symbol_count(seq, 'N') AS n
    FROM butterflies.regions
)
SELECT *, round(100.0 * (g + c) / NULLIF(a + c + g + t, 0), 4)
          AS gc_percent
FROM counts
ORDER BY key;
