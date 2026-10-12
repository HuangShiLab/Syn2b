#!/bin/bash
#SBATCH -J hrom-s2b-frag
#SBATCH -p intel
#SBATCH -o /lustre1/g/aos_shihuang/syn2b_benchmarks/logs/hrom_fragmentation_%j.out
#SBATCH -e /lustre1/g/aos_shihuang/syn2b_benchmarks/logs/hrom_fragmentation_%j.err
#SBATCH -c 16
#SBATCH --mem=140G
#SBATCH -t 24:00:00
set -euo pipefail
ROOT=/lustre1/g/aos_shihuang/syn2b_benchmarks
INPUT=$ROOT/inputs/hrom
OUT=$ROOT/results/hrom_full_fragmentation
SYN2B=$ROOT/Syn2b/target/release/syn2b
mkdir -p "$OUT"
/usr/bin/time -v python3 "$ROOT/Syn2b/scripts/benchmark_contig_fragmentation.py" \
  --reference-tgt-dir /lustre1/g/aos_shihuang/databases/HROM/census/tgt \
  --metadata "$INPUT/hrom_metadata.tsv" \
  --query-list "$INPUT/query_list.tsv" \
  --syn2b "$SYN2B" \
  --output "$OUT/fragment_hits.tsv" \
  --summary-output "$OUT/fragment_summary.tsv" \
  --workdir "$OUT/work" \
  --mode 2brad --enzymes BcgI \
  --min-tags 5 --min-adjacencies 3 --max-targets 25 \
  2>&1 | tee "$OUT/runtime.txt"
