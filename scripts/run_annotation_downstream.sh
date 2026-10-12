#!/usr/bin/env bash
set -euo pipefail
root=/Users/macstudio/Downloads/Syn2b
func=$root/external_data/oapgc_functional
syn2b=$root/target/release/syn2b
# Wait for eggNOG output and DIAMOND to exit.
while [[ ! -s $func/oapgc_representatives.emapper.annotations ]] || pgrep -f 'diamond blastp.*oapgc_representative_proteins' >/dev/null; do
  sleep 300
done
# Wait for antiSMASH queue to exit.
while pgrep -f "$root/scripts/run_antismash_representatives.sh" >/dev/null || pgrep -f 'antismash/standalone:8.0.4' >/dev/null; do
  sleep 300
done
# Retry any failed antiSMASH jobs once.
find $root/external_data/oapgc_antismash -name '.failed' -delete
bash $root/scripts/run_antismash_representatives.sh
# Build feature tables and searchable indexes.
bash $root/scripts/build_oapgc_functional_clusters.sh
# MIBiG retrieval: OAPGC antiSMASH clusters query MIBiG index.
$syn2b cluster-build --features $func/oapgc_antismash_features.tsv \
  --output $func/oapgc_antismash_query.s2cl --min-genes 1
$syn2b search-cluster --query $func/oapgc_antismash_features.tsv \
  --index $root/external_data/imported/clusters/mibig.s2cl \
  --output $func/mibig_retrieval.tsv \
  --min-functional-jaccard 0.10 --min-synteny 0.05 --max-targets 10
# CARD query against AMRFinderPlus + ResFinder references.
$syn2b cluster-build --features $root/external_data/imported/card_features.tsv \
  --output $root/external_data/imported/clusters/card_query.s2cl --min-genes 1
$syn2b cluster-index --cluster-files \
  $root/external_data/imported/clusters/amrfinder.s2cl \
  $root/external_data/imported/clusters/resfinder.s2cl \
  --output $root/external_data/imported/amrfinder_resfinder.s2cldb
$syn2b search-cluster --query $root/external_data/imported/card_features.tsv \
  --index $root/external_data/imported/amrfinder_resfinder.s2cldb \
  --output $func/card_amr_crosswalk.tsv \
  --min-functional-jaccard 0.10 --min-synteny 0.0 --max-targets 10
echo "$(date '+%F %T') downstream annotation pipeline complete" > $func/downstream.complete
