#!/usr/bin/env bash
set -u
input=$1
root=/Users/macstudio/Downloads/Syn2b
gid=$(basename "$input" .fna)
outdir=$root/external_data/oapgc_antismash/$gid
logdir=$root/external_data/oapgc_antismash/logs
mkdir -p "$logdir"
marker=$outdir/.complete
if [[ -e $marker ]]; then exit 0; fi
# A failed run may leave partial files; archive it so antiSMASH gets a clean dir.
if [[ -d $outdir ]]; then
  archive=$root/external_data/oapgc_antismash/failed_runs/$gid-$(date +%Y%m%d-%H%M%S)
  mkdir -p "$(dirname "$archive")"
  mv "$outdir" "$archive"
fi
mkdir -p "$outdir"
docker run --rm --platform linux/amd64 \
  -v $root/external_data/oapgc_clean_fasta_representatives_fna:/input:ro \
  -v $root/external_data/oapgc_antismash:/output \
  --entrypoint /usr/local/bin/antismash \
  antismash/standalone:8.0.4 \
  /input/$(basename "$input") \
  --output-dir /output/$gid \
  --output-basename $gid \
  --taxon bacteria \
  --minimal \
  --genefinding-tool prodigal \
  --cb-knownclusters \
  --cpus 2 \
  > $logdir/$gid.log 2>&1
status=$?
if [[ $status -eq 0 ]]; then touch $marker; else echo $status > $outdir/.failed; fi
exit $status
