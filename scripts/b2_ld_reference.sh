#!/usr/bin/env bash
# 1000G Europeans, four gene windows: counts, freqs
set -euo pipefail
D=~/o_data
W=$D/work
PROJ=${PROJ:-$(cd "$(dirname "$0")/.." && pwd)}
P2="$D/bin/plink2 --memory 6500 --threads 16"
OUT=$PROJ/data/ld
mkdir -p $OUT $W
awk -F'\t' 'NR>1 && $5=="EUR"{print $1}' $D/all_phase3.psam > $W/eur.txt
wc -l < $W/eur.txt
tail -n +2 $PROJ/inputs/genes.csv | tr -d '\r' | while IFS=, read -r gene chr s e rest; do
  $P2 --pfile $D/all_phase3 vzs --keep $W/eur.txt --chr $chr \
    --from-bp $((s - 100000)) --to-bp $((e + 100000)) --snps-only just-acgt --max-alleles 2 \
    --maf 0.01 --set-all-var-ids '@:#:$r:$a' --rm-dup exclude-all \
    --export A --freq --out $OUT/$gene > /dev/null
  echo "$gene $(head -1 $OUT/$gene.raw | wc -w)"
done
