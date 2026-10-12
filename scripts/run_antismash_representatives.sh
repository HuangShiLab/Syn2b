#!/usr/bin/env bash
set -u
root=/Users/macstudio/Downloads/Syn2b
find $root/external_data/oapgc_clean_fasta_representatives_fna -maxdepth 1 -type f -name '*.fna' -print0 \
  | xargs -0 -P ${ANTISMASH_JOBS:-2} -n 1 $root/scripts/run_one_antismash.sh
