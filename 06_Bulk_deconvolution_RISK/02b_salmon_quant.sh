#!/bin/bash
# Per sample: download FASTQ (md5-checked), salmon quant, delete FASTQ
# Usage: 02b_salmon_quant.sh <index dir> <out dir> <ncores> <rows of risk_runs.tsv>...
set -euo pipefail
idx=$1; out=$2; ncores=$3; shift 3
for row in "$@"; do
  IFS=$'\t' read -r run gsm sample_id layout reads md5 ftp <<< "$row"
  [[ -s "$out/$sample_id/quant.sf" ]] && continue   # already quantified
  for try in 1 2 3 4 5; do
    wget -nv --tries=3 --waitretry=30 -O "$run.fastq.gz" "https://$ftp" &&
      echo "$md5  $run.fastq.gz" | md5sum -c --quiet && break
    echo "$run: download attempt $try failed"; rm -f "$run.fastq.gz"; sleep 60
  done
  [[ -s "$run.fastq.gz" ]] || { echo "$run: download failed"; exit 1; }
  salmon quant -i "$idx" -l A -r "$run.fastq.gz" -p "$ncores" -o "$sample_id"
  rm -f "$run.fastq.gz"
  mkdir -p "$out" && cp -r "$sample_id" "$out/" && rm -rf "$sample_id"
done
