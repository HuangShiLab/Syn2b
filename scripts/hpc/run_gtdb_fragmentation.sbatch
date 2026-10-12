#!/bin/bash
#SBATCH -J gtdb-s2b-frag
#SBATCH -p intel
#SBATCH -o /lustre1/g/aos_shihuang/syn2b_benchmarks/logs/gtdb_fragmentation_%j.out
#SBATCH -e /lustre1/g/aos_shihuang/syn2b_benchmarks/logs/gtdb_fragmentation_%j.err
#SBATCH -c 16
#SBATCH --mem=180G
#SBATCH -t 48:00:00
set -euo pipefail
ROOT=/lustre1/g/aos_shihuang/syn2b_benchmarks
INPUT=$ROOT/inputs/gtdb
OUT=$ROOT/results/gtdb_full_fragmentation
SYN2B=$ROOT/Syn2b/target/release/syn2b
mkdir -p "$OUT"
/usr/bin/time -v python3 "$ROOT/Syn2b/scripts/benchmark_contig_fragmentation.py" \
  --reference-tgt-dir "$INPUT/tgt" \
  --metadata "$INPUT/gtdb_metadata.tsv" \
  --query-list "$INPUT/query_list.tsv" \
  --syn2b "$SYN2B" \
  --output "$OUT/fragment_hits.tsv" \
  --summary-output "$OUT/fragment_summary.tsv" \
  --workdir "$OUT/work" \
  --mode 2brad --enzymes BcgI \
  --min-tags 5 --min-adjacencies 3 --max-targets 25 \
  2>&1 | tee "$OUT/runtime.txt"
