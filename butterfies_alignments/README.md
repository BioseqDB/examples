# Butterfly alignments

Code for [Comparing butterfly genomes with BioseqDB, SQLAlchemy and Plotly](https://bioseqdb.com/blog/comparing-butterfly-genomes-with-bioseqdb/).

Sequences are retrieved from NCBI by BioseqDB SQL, stored locally, and aligned inside PostgreSQL. Python submits those queries, checks their results and draws figures. There is no Minimap2 or R dependency. Downloaded sequences and annotations are not committed.

## Database setup

**Availability, 16 October 2026:** external-database support is not enabled in the published `0.1.0-alpha.1` package build profiles. The PostgreSQL 17 ARM64 Docker image, Ubuntu package and macOS package were checked directly. This is planned for the next alpha release. The tutorial currently needs an external-enabled source build; it cannot be enabled in those packages with a setting or by creating the missing tables.

The [usual installation instructions](https://bioseqdb.com/v0.1/tutorials/getting_started/) cover the released packages. Before using one for this tutorial, connect as its database owner and check:

```sql
SELECT to_regnamespace('bioseq_external') IS NOT NULL AS external_schema,
       to_regclass('bioseq_external.catalog') IS NOT NULL AS catalog,
       to_regclass('bioseq_external.sequences') IS NOT NULL AS sequences;
```

All three must be true. If the feature is present but your role cannot read its tables, the database owner can grant `USAGE ON SCHEMA bioseq_external` and `SELECT ON bioseq_external.catalog, bioseq_external.sequences` to that role. It also needs permission to create the `butterflies` schema and its own tables. Do not grant superuser merely for this example.

The source-build option is `-DBIOSEQ_BUILD_EXTERNAL_API_CLIENTS=ON`. The provided `build_macos.sh` documents the tested Apple Silicon / PostgreSQL 17 / LLVM 22 configuration. It expects a BioseqDB source checkout and a bootstrapped vcpkg installation. The small `ncbi-provider.patch` adds organism-name catalog selection and bounded NCBI sequence fetching to the pinned public source revision; skip that patch only if your checkout already contains those changes. See the script for the exact commands and prerequisites.

From this example directory, prepare that pinned source in a separate checkout:

```sh
export BIOSEQ_SOURCE="$PWD/../../bioseqdb-tutorial-source"
git clone --recurse-submodules https://github.com/styczynski/bioseqdb.git "$BIOSEQ_SOURCE"
git -C "$BIOSEQ_SOURCE" checkout f7a34a198e42ffa345917eac9f02e41116d34fe7
git -C "$BIOSEQ_SOURCE" submodule update --init --recursive
git -C "$BIOSEQ_SOURCE" apply "$PWD/ncbi-provider.patch"
bash build_macos.sh
cmake --install "$BIOSEQ_SOURCE/build/butterflies-tutorial" --component bioseq
```

Export `VCPKG_ROOT` first, as described at the top of `build_macos.sh`. Building the C++ module graph takes considerably longer than the analysis itself. The helper does not start, stop or overwrite a running PostgreSQL server. The final install command targets the PostgreSQL selected by the build, so use a dedicated development installation.

Once the source build is installed, start a separate development cluster. These commands use a new directory and port 55432, and leave any existing server running. `initdb` asks you to choose a password for its `postgres` role:

```sh
BIOSEQ_TUTORIAL_PG="$(brew --prefix postgresql@17)"
BIOSEQ_TUTORIAL_DATA="$PWD/.postgres-butterflies"
"$BIOSEQ_TUTORIAL_PG/bin/initdb" -D "$BIOSEQ_TUTORIAL_DATA" \
  -U postgres --auth=scram-sha-256 --pwprompt
cat >> "$BIOSEQ_TUTORIAL_DATA/postgresql.conf" <<'CONFIG'
listen_addresses = 'localhost'
port = 55432
shared_preload_libraries = 'bioseq_cache,bioseq'
CONFIG
"$BIOSEQ_TUTORIAL_PG/bin/pg_ctl" -D "$BIOSEQ_TUTORIAL_DATA" \
  -l "$BIOSEQ_TUTORIAL_DATA/server.log" start
"$BIOSEQ_TUTORIAL_PG/bin/createdb" -h localhost -p 55432 -U postgres butterflies
"$BIOSEQ_TUTORIAL_PG/bin/psql" -h localhost -p 55432 -U postgres -d butterflies \
  -v ON_ERROR_STOP=1 -c 'CREATE EXTENSION bioseq_cache; CREATE EXTENSION bioseq;'
```

The server must be able to reach NCBI over HTTPS. Replace the demonstration password in `DATABASE_URL` below with the one you chose. When finished, stop only this cluster with `"$BIOSEQ_TUTORIAL_PG/bin/pg_ctl" -D "$BIOSEQ_TUTORIAL_DATA" stop`.

## Run

Use Python 3.11–3.13 and Poetry 2. If you already cloned the examples for the source-build steps, skip the first two commands:

```sh
git clone https://github.com/BioseqDB/examples.git
cd examples/butterfies_alignments
poetry install
export DATABASE_URL='postgresql+psycopg://postgres:bioseq@localhost:55432/butterflies'
poetry run butterfly-tutorial
```

The command creates `butterflies.regions`, replaces only the two tutorial keys, and writes generated results under `output/`. The imports are transactional. On a later run, `poetry run butterfly-tutorial --reuse` skips the NCBI sequence imports but verifies their hashes again. Annotation files have a checked local cache. Use `--no-figures` to run SQL and validation without static image export.

Plotly's SVG export uses Kaleido, which needs Chrome. If Chrome is absent, install it normally or run `poetry run plotly_get_chrome`. The interactive HTML files share a local `plotly.min.js`; keep that file beside them when copying the figures elsewhere. Each HTML figure has a light/dark toggle and also follows the initial browser colour preference.

```sh
poetry run jupyter lab notebooks/butterfly_alignments.ipynb
```

The notebook contains the tutorial's Python snippets and runs the same SQL files. Execute it from the Poetry environment after exporting `DATABASE_URL`. Its committed copy has no saved outputs or credentials.

## Inputs and checks

`sources.json` pins the sequence accession versions, one-million-base intervals, GFF URLs, SHA-256 checksums and annotation features. All genomic coordinates there are zero-based and half-open. The erato record is an **unplaced scaffold**, and its submitted annotation calls the optix homolog `AC2CFC_001707` / homeobox protein SIX6; that naming distinction is retained in the plots and source metadata.

The runner checks every imported length, sequence hash, alphabet, base-count sum, alignment coordinate and CIGAR operation. It replays all alignment blocks against the actual bases and checks the reported edit counts. It also extracts the annotated coding sequences in transcript orientation and checks their nucleotide and amino-acid comparison. Empty imports, empty alignments and empty plotting selections fail explicitly.

The plotting filter retains at least 300 alignment columns, mapping quality at least 20 and identity at least 70%. It does not alter the saved unfiltered matches. Identity counts gap columns as differences. Genomic span coverage merges overlapping intervals; it is not the sum of alignment lengths.

Gene annotations and measured GC tracks belong to the current NCBI assemblies. No historical TE or ATAC tracks are relabelled as measurements on these references. The ATAC diagram is explanatory only.

## Files

- `sql/`: discovery, schema, NCBI imports, inspection, composition, alignment and windowed GC queries.
- `butterfly_tutorial/workflow.py`: SQLAlchemy execution and independent validation.
- `butterfly_tutorial/annotations.py`: checked NCBI GFF download and feature extraction.
- `butterfly_tutorial/plots.py`: empirical Plotly figures and HTML/SVG export.
- `butterfly_tutorial/diagrams.py`: original explanatory SVG illustrations.
- `figures/`: the article's measured charts and explanatory diagrams as SVGs, with light and dark palettes.
- `notebooks/butterfly_alignments.ipynb`: executable notebook.
- `schema.tex`: LaTeX counterpart of the table field descriptions.
- `validation/`: small measured summaries and review records, with no input sequences.

The long-form explanation stays on the website. Follow the accession and source links in `sources.json` for the biological data provenance.
