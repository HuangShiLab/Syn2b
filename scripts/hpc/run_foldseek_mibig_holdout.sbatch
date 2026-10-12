#!/bin/bash
#SBATCH -J s2b-fold-holdout
#SBATCH -p intel
#SBATCH -o /lustre1/g/aos_shihuang/syn2b_benchmarks/logs/foldseek_holdout_%j.out
#SBATCH -e /lustre1/g/aos_shihuang/syn2b_benchmarks/logs/foldseek_holdout_%j.err
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
set -euo pipefail
ROOT=/lustre1/g/aos_shihuang/syn2b_benchmarks
FOLDSEEK=/lustre1/g/aos_shihuang/tools/foldseek/bin/foldseek
OUT=$ROOT/results/foldseek_mibig_holdout
mkdir -p "$OUT/tmp"
/usr/bin/time -v "$FOLDSEEK" easy-search \
  "$ROOT/input_proteins/heldout_proteins.faa" \
  "$ROOT/input_proteins/target_proteins.faa" \
  "$OUT/heldout_vs_target.fold.m8" "$OUT/tmp" \
  --threads 16 --alignment-type 1 \
  --format-output query,target,pident,alnlen,evalue,bits,qstart,qend,tstart,tend,qcov,tcov \
  2>&1 | tee "$OUT/runtime.txt"
