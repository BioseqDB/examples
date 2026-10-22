# Scoring matrices from 1,000 human mitochondrial genomes

Code for [Scoring matrices from human genomes](https://bioseqdb.com/blog/scoring-matrices-from-human-genomes/).

This example downloads 1,000 mitochondrial genomes from NCBI through BioseqDB, aligns them to the revised Cambridge Reference Sequence (rCRS), and counts aligned nucleotide and amino-acid pairs in PostgreSQL. SQL also calculates log-odds scores, substitution spectra and distances. Python reads the exported tables and draws the figures.

The accession manifest selects the first 1,000 records titled `complete genome` from the 2,000-person Danish study by [Li et al. (2014)](https://doi.org/10.1038/ejhg.2013.282). Each record has a different submitted isolate identifier. The reproducible accession subset is described in `samples.csv` and `sources.json`. These are complete mitochondrial genome records, roughly 16.6 kb per person, with some uncalled bases represented by `N`.

## Database setup

Use PostgreSQL with BioseqDB built with external-database support. The source-build procedure in [the butterfly example](https://github.com/BioseqDB/examples/tree/main/butterfies_alignments#database-setup) describes the tested configuration.

**Availability on 16 October 2026:** the published `0.1.0-alpha.1` Docker images and native packages omit the external-database integration. This feature is planned for the next alpha release. A source build must enable `BIOSEQ_BUILD_EXTERNAL_API_CLIENTS=ON`; a runtime setting cannot add code missing from the binary.

The database needs the `bioseq` and `bioseq_cache` extensions. Use a database owner to install them, or ask the owner to do so before running the example. Your role needs permission to create the `scoring` schema and read `bioseq_external.sequences`. Confirm the external table exists:

```sql
SELECT to_regclass('bioseq_external.sequences') IS NOT NULL
       AS external_sequences_available;
```

## Run the analysis

The database must already exist. For a new source-build server on port 55432, create it once with `createdb -h localhost -p 55432 -U postgres scoring_matrices`. Run this example in its own database or use an empty `scoring` schema. The PostgreSQL client tools `psql` and `createdb` must be on your `PATH`.

```sh
git clone https://github.com/BioseqDB/examples.git
cd examples/scoring_matrices
export DATABASE_URL='postgresql://postgres:your-password@localhost:55432/scoring_matrices'
bash run.sh
```

`DATABASE_URL` is a PostgreSQL URI or connection string accepted by `psql`. Keep your password outside version control.

If you downloaded the ZIP, extract it and run the commands from its `scoring_matrices` directory, starting with the `DATABASE_URL` export.

Downloads are sequential, with a delay of at least 0.4 seconds between NCBI requests. The import skips accessions already in the table, so rerunning `bash run.sh` resumes after an interrupted download. Reruns reload the accession metadata and replace this example's derived tables in the `scoring` schema. Downloaded genome sequences remain in your database and are excluded from this repository.

The SQL runner writes small result tables to `results/`. Use Python 3.11 or later to draw the figures:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python render_charts.py
```

The plotting dependency is Matplotlib 3.10.7. SQL performs the alignment, sequence translation and matrix estimation. The Python script only plots its exported results.

## What the counts mean

Every sample is aligned to `NC_012920.1`. CIGAR operations place each called base at a reference position, so insertions and deletions cannot shift all later comparisons. Insertions are counted separately. Deletions and positions outside an alignment are represented by `-` in the projected sequence. In the measured run, 999 records cover the entire reference. `KF161455.1` covers every query base but leaves seven terminal reference positions uncovered; those positions are masked.

`N` means an uncalled nucleotide and is excluded from pair counts. The rCRS also contains a historical `N` spacer at position 3107; the whole reference position is excluded from nucleotide substitution calculations. Raw `N` counts and aligned gaps remain in the quality reports.

The directional reference spectrum counts differences from the rCRS letter to the sample letter. Its arrows describe this comparison direction. The symmetric cohort matrix compares all different samples at each aligned position. It uses per-position letter counts, which gives the same counts as enumerating every pair without performing 499,500 pairwise alignments.

The protein comparison uses all 13 rCRS coding sequences and NCBI mitochondrial genetic code 2. SQL reverse-complements ND6 and applies the mitochondrial initiation-codon rule. An entire sample/CDS is excluded if it contains an insertion, deletion, incomplete alignment, premature stop, invalid initiation codon or ambiguous initiation codon. Within retained CDSs, codons containing `N` are excluded. Protein `N` means asparagine and remains part of the 20-letter alphabet. Partial terminal stop codons are outside the reference protein lengths.

The empirical matrices use observed frequencies without pseudocounts. A pair with no observations has a SQL `NULL` log-odds score and is shown as unobserved in the figures. The conditional difference spectrum omits diagonal matches before normalizing each row. BLOSUM62 and PAM250 probabilities come from NCBI's public-domain reference data, retained in `reference_matrices/` with their source notice.

## Files

- `sql/01_schema.sql`: creates the sequence table and loads sample metadata.
- `sql/02_import.sql`: downloads the versioned NCBI records through BioseqDB.
- `sql/00_preflight.sql`: checks that the complete cohort imported before analysis.
- `sql/03_quality.sql`: reports lengths and unknown bases.
- `sql/04_align.sql`: performs BioseqDB alignments and stores CIGARs.
- `sql/05_columns.sql`: projects the alignment and counts each letter per site.
- `sql/06_nucleotide_matrix.sql`: estimates nucleotide pair scores and reference differences.
- `sql/07_proteins.sql`: translates coding sequences and counts amino-acid pairs.
- `sql/08_reference_matrices.sql`: loads published matrix probabilities and calculates comparable spectra.
- `sql/10_distance.sql`: calculates pairwise p-distances for eight selected samples.
- `sql/11_validate.sql`: checks the full cohort, CIGAR replay, reference translation and matrix arithmetic.
- `sql/09_export.sql`: exports the plotted result tables.

`sql/00_toy_examples.sql` contains the introductory log-odds and column-count examples from the article. It can run independently with `psql "$DATABASE_URL" -X --set ON_ERROR_STOP=on --file sql/00_toy_examples.sql`. The committed `results/` files contain the measured summary tables, so `python render_charts.py` can also redraw the figures before downloading the cohort.

The validation step requires exactly 1,000 unique samples, full query coverage and at least 99.9% reference coverage, agreement with every reference protein annotation, matching CIGAR edit counts and normalized symmetric pair matrices. It stops the run when any of these checks fails.
