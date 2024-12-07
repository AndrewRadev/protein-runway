#! /bin/sh

if [ $# -lt 1 ]; then
  echo "USAGE: bash scripts/generate_snakemake_graphs.sh <protein-name>"
  exit 1
fi

protein_name=$1

mkdir -p 04_extra/graphs/

snakemake --dag 03_output/$protein_name.segmentation.tsv \
  | dot -Tsvg > 04_extra/graphs/dag.svg
snakemake --filegraph 03_output/$protein_name.segmentation.tsv \
  | dot -Tsvg > 04_extra/graphs/filegraph.svg
