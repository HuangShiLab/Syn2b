#!/usr/bin/env bash
set -euo pipefail
root=/Users/macstudio/Downloads/Syn2b
func=$root/external_data/oapgc_functional
syn2b=$root/target/release/syn2b
python3 $root/scripts/emapper_to_features.py \
  --annotation-dir /Users/macstudio/Downloads/OAPGC/annotations_raw \
  --emapper $func/oapgc_representatives.emapper.annotations \
  --output $func/oapgc_eggnog_features.tsv
python3 $root/scripts/antismash_to_features.py \
  --antismash-dir $root/external_data/oapgc_antismash \
  --output $func/oapgc_antismash_features.tsv
$syn2b cluster-build --features $func/oapgc_antismash_features.tsv \
  --output $func/oapgc_antismash.s2cl --min-genes 1
$syn2b cluster-index --cluster-files $func/oapgc_antismash.s2cl \
  --output $func/oapgc_antismash.s2cldb
# Keep the complete gene-level annotation table as a separate search artifact.
$syn2b cluster-build --features $func/oapgc_eggnog_features.tsv \
  --output $func/oapgc_eggnog_genes.s2cl --min-genes 1
$syn2b cluster-index --cluster-files $func/oapgc_eggnog_genes.s2cl \
  --output $func/oapgc_eggnog_genes.s2cldb
ls -lh $func/*features.tsv $func/*.s2cl $func/*.s2cldb
