SELECT t.accession AS target_accession,
       q.accession AS query_accession,
       bioseq.target_start(t.seq) AS target_start,
       bioseq.target_end(t.seq) AS target_end,
       bioseq.query_start(t.seq) AS query_start,
       bioseq.query_end(t.seq) AS query_end,
       bioseq.is_reverse(t.seq) AS is_reverse,
       bioseq.score(t.seq) AS score,
       bioseq.cigar(t.seq) AS cigar,
       bioseq.edit(t.seq) AS edits,
       bioseq.block_length(t.seq) AS block_length,
       bioseq.mapq(t.seq) AS mapq,
       bioseq.is_primary(t.seq) AS is_primary
FROM butterflies.regions AS t
JOIN butterflies.regions AS q
  ON mmsearch(t.seq, q.seq, aligner_options(
      kmer_size => 15, window_size => 10,
      max_gap => 200, max_gap_ref => 200, min_chain_score => 40))
WHERE t.key = 'melpomene' AND q.key = 'erato'
ORDER BY score DESC, target_start, query_start;
