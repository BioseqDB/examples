WITH windows AS (
    SELECT r.accession, w.start,
           least(w.start + :window_size, bioseq.len(r.seq)) AS "end",
           r.seq[w.start:least(w.start + :window_size, bioseq.len(r.seq))]
               AS piece
    FROM butterflies.regions AS r
    CROSS JOIN LATERAL generate_series(
        0, bioseq.len(r.seq) - 1, :window_size) AS w(start)
), counts AS (
    SELECT accession, start, "end", bioseq.len(piece) AS length,
           bioseq.symbol_count(piece, 'G') + bioseq.symbol_count(piece, 'C') AS gc,
           bioseq.symbol_count(piece, 'N') AS n
    FROM windows
)
SELECT accession, start, "end",
       round(100.0 * gc / NULLIF(length - n, 0), 4) AS gc_percent,
       1.0 * n / length AS n_fraction
FROM counts
ORDER BY accession, start;
