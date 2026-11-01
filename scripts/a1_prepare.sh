#!/usr/bin/env bash
# Genotype QC, pruning, PCs, kinship, pool, test set
# WSL, big files in ~/o_data. MSI: BASE_ONLY=1
# rebuilds only the QC'd base (rest is shipped)
set -euo pipefail
PROJ=${PROJ:-$(cd "$(dirname "$0")/.." && pwd)}
D=${O_DATA:-$HOME/o_data}
P2="$D/bin/plink2 --memory ${MEM:-6500} --threads ${THREADS:-16}"
W=$D/work
mkdir -p $W
cd $D

if [ ! -f all_phase3.pgen ]; then
  $D/bin/plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
fi
[ -f all_phase3.psam ] || cp phase3_corrected.psam all_phase3.psam
[ -f all_phase3.pvar.zst ] || cp all_phase3_noannot.pvar.zst all_phase3.pvar.zst

# biallelic autosomal SNPs, no missing, MAF >= 0.01
# per chromosome: 85M-variant pvar too big for WSL
rm -f $W/merge.txt
for c in $(seq 1 22); do
  $P2 --pfile all_phase3 vzs --chr $c --snps-only just-acgt --max-alleles 2 \
    --geno 0 --maf 0.01 --set-all-var-ids '@:#:$r:$a' --rm-dup exclude-all \
    --make-pgen --out $W/base_chr$c > /dev/null
  echo $W/base_chr$c >> $W/merge.txt
done
$P2 --pmerge-list $W/merge.txt pfile --make-pgen --out $W/base > /dev/null
rm -f $W/base_chr* $W/base-merge.*
$P2 --pfile $W/base --freq --out $W/base > /dev/null
[ "${BASE_ONLY:-0}" = 1 ] && exit 0

$P2 --pfile $W/base --maf 0.05 --indep-pairwise 1000kb 1 0.2 --out $W/prune
$P2 --pfile $W/base --extract $W/prune.prune.in --pca 10 --out $W/pcs
$P2 --pfile $W/base --extract $W/prune.prune.in --make-king-table \
  --king-table-filter 0.0884 --out $W/king

# pool: chr 1-11, MAF >= 0.05 (filter before thin)
$P2 --pfile $W/base --chr 1-11 --maf 0.05 --write-snplist --out $W/pool_all
$P2 --pfile $W/base --extract $W/pool_all.snplist --thin-count 20000 --seed 2027 \
  --export A --out $W/pool
# test set: chr 12-22, all MAF >= 0.01 SNPs, thinned
$P2 --pfile $W/base --chr 12-22 --thin-count 300000 --seed 2027 \
  --write-snplist --out $W/testset

mkdir -p $PROJ/data/work
cp $W/pool.raw $D/all_phase3.psam $W/pcs.eigenvec $W/pcs.eigenval $W/king.kin0 $PROJ/data/work/
wc -l $W/base.pvar $W/prune.prune.in $W/testset.snplist $W/king.kin0
