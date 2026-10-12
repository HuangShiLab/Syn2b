#!/bin/bash
#SBATCH -J gtdb-s2b-digest
#SBATCH -p intel
#SBATCH -o /lustre1/g/aos_shihuang/syn2b_benchmarks/logs/gtdb_digest_%A_%a.out
#SBATCH -e /lustre1/g/aos_shihuang/syn2b_benchmarks/logs/gtdb_digest_%A_%a.err
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -a 0-28%29
set -euo pipefail
ROOT=/lustre1/g/aos_shihuang/syn2b_benchmarks
LIST=/home/shihuang/fast2brad_holo/gtdb_representative_paths.tsv
OUT=$ROOT/gtdb/tgt
SYN2B=$ROOT/Syn2b/target/release/syn2b
CHUNK=5000
START=$((SLURM_ARRAY_TASK_ID * CHUNK + 1))
END=$(((SLURM_ARRAY_TASK_ID + 1) * CHUNK))
mkdir -p "$OUT" "$ROOT/logs"
tail -n +$START "$LIST" | head -n $CHUNK | while IFS=$'\t' read -r species path; do
  base=$(basename "$path"); base=${base%.fna.gz}; base=${base%.fna}
  dst="$OUT/$base.tgt"
  if [[ ! -s "$dst" ]]; then
    "$SYN2B" digest --input "$path" --output "$dst" --format text --mode 2brad --enzymes BcgI
  fi
done
