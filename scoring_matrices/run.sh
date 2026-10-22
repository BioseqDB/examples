#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results
: "${DATABASE_URL:?Set DATABASE_URL to a PostgreSQL connection URI or conninfo string}"
for step in 01_schema 02_import 00_preflight 03_quality 04_align 05_columns 06_nucleotide_matrix 07_proteins 08_reference_matrices 10_distance 11_validate 09_export; do
  printf '\nRunning %s\n' "$step"
  psql "$DATABASE_URL" -X --set ON_ERROR_STOP=on --file "sql/$step.sql"
done
psql "$DATABASE_URL" -X --set ON_ERROR_STOP=on --file tests/invariants.sql
