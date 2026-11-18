#!/bin/bash
# Salmon index: GENCODE v44 transcripts, GRCh38 primary assembly as decoys; tx2gene table
# Usage: 02a_salmon_index.sh <index dir> <ncores>
set -euo pipefail
out=$1; ncores=$2
base=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44
wget -q "$base/gencode.v44.transcripts.fa.gz" "$base/GRCh38.primary_assembly.genome.fa.gz"
zcat GRCh38.primary_assembly.genome.fa.gz | grep '^>' | cut -d ' ' -f 1 | sed 's/^>//' > decoys.txt
cat gencode.v44.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz
zcat gencode.v44.transcripts.fa.gz | awk -F'|' '/^>/{sub(/^>/, "", $1); print $1 "\t" $2 "\t" $6}' > tx2gene.tsv
salmon index -t gentrome.fa.gz -d decoys.txt -i "$out" --gencode -p "$ncores"
cp tx2gene.tsv "$out/"
