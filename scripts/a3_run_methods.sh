#!/usr/bin/env bash
# M1, M2, M2f, M3, M3s (+ K3), condensed per batch
# PART=setup|q|b|all; TRAITS=a,b overrides traits
# SMOKE=1: two traits per kind, for timing
# laptop defaults; msi/ sets O_DATA, INPUTS, bins
set -euo pipefail
PROJ=${PROJ:-$(cd "$(dirname "$0")/.." && pwd)}
O_DATA=${O_DATA:-$HOME/o_data}
INPUTS=${INPUTS:-$PROJ/data/work}
THREADS=${THREADS:-16}
W=$O_DATA/work
P2="${PLINK2_BIN:-$O_DATA/bin/plink2} --memory ${MEM:-6500} --threads $THREADS"
# libgomp from the Ubuntu .deb, no root (laptop)
[ -d $O_DATA/lib/x ] && export LD_LIBRARY_PATH=$O_DATA/lib/x/usr/lib/x86_64-linux-gnu
RG="${REGENIE_BIN:-$O_DATA/bin/mkl/regenie_v4.1.3.gz_x86_64_Linux_mkl} --threads $THREADS"
export O_WORK=$W
SUM="python3 $PROJ/scripts/a4_summarize.py"
PART=${PART:-all}
TAG=${TAG:-main}
O=$W/out_$TAG
mkdir -p $O
cd $W

setup() {
  cp $INPUTS/pheno_q.txt $INPUTS/pheno_b.txt $INPUTS/truth.csv $INPUTS/large_snps.txt .
  cut -f2- pheno_q.txt | sed '1s/^IID/#IID/' > pheno_q_iid.txt
  cut -f2- pheno_b.txt | sed '1s/^IID/#IID/' > pheno_b_iid.txt
  awk 'NR==1{printf "FID IID"; for(i=2;i<=NF;i++) printf " %s",$i; print ""; next}
       {printf "%s %s",$1,$1; for(i=2;i<=NF;i++) printf " %s",$i; print ""}' pcs.eigenvec > pcs_fidiid.txt
  cat testset.snplist large_snps.txt | sort -u > tests.snplist
  # REGENIE wants #FID in .psam; link the genotypes
  ln -sf base.pgen rg.pgen
  ln -sf base.pvar rg.pvar
  awk 'NR==1{print "#FID\tIID" substr($0, index($0, "\t")); next}{print $1 "\t" $0}' base.psam > rg.psam
  python3 -c "
import random; random.seed(2027)
ids = open('testset.snplist').read().split()
open('keep200.txt','w').write('\n'.join(random.sample(ids, 200)) + '\n')"
  $P2 --pfile base --extract keep200.txt --export A --out keep200 > /dev/null
}

if [ "${SMOKE:-0}" = 1 ]; then
  QN=Q01,Q02; BN=B20_01,B02_01; RES=${RES:-$INPUTS/smoke.jsonl}
else
  QN=$(head -1 $INPUTS/pheno_q.txt | cut -f3- | tr '\t' ','); BN=$(head -1 $INPUTS/pheno_b.txt | cut -f3- | tr '\t' ',')
  RES=${RES:-$INPUTS/summary.jsonl}
fi
if [ -n "${TRAITS:-}" ]; then QN=$TRAITS; BN=$TRAITS; fi

stamp() { echo "[$(date +%T)] $1 (${SECONDS}s)"; }
condense() { # method format files...
  local m=$1 f=$2; shift 2
  $SUM $m $f $RES "$@" > /dev/null
  rm -f "$@"
  stamp "$m done"
}

run_q() {
  $P2 --pfile base --extract tests.snplist --pheno pheno_q_iid.txt --pheno-name $QN \
    --glm allow-no-covars hide-covar --out $O/m1q > /dev/null
  condense M1_q plink $O/m1q.*.glm.linear
  $P2 --pfile base --extract tests.snplist --pheno pheno_q_iid.txt --pheno-name $QN \
    --covar pcs.eigenvec --covar-variance-standardize --glm hide-covar --out $O/m2q > /dev/null
  condense M2_q plink $O/m2q.*.glm.linear
  $RG --step 1 --pgen rg --extract prune.prune.in --phenoFile pheno_q.txt --phenoColList $QN \
    --covarFile pcs_fidiid.txt --bsize 1000 --lowmem --lowmem-prefix $O/tmpq --out $O/s1q > /dev/null
  stamp "REGENIE step 1 quantitative"
  $RG --step 2 --pgen rg --extract tests.snplist --phenoFile pheno_q.txt --phenoColList $QN \
    --covarFile pcs_fidiid.txt --pred $O/s1q_pred.list --bsize 400 --out $O/m3q > /dev/null
  condense M3_q regenie $O/m3q_*.regenie
  $RG --step 2 --pgen rg --extract tests.snplist --phenoFile pheno_q.txt --phenoColList $QN \
    --covarFile pcs_fidiid.txt --ignore-pred --bsize 400 --out $O/k3q > /dev/null
  condense K3_q regenie $O/k3q_*.regenie
}

run_b() {
  $P2 --pfile base --extract tests.snplist --pheno pheno_b_iid.txt --pheno-name $BN --1 \
    --glm allow-no-covars hide-covar no-firth --out $O/m1b > /dev/null
  condense M1_b plink $O/m1b.*.glm.logistic
  $P2 --pfile base --extract tests.snplist --pheno pheno_b_iid.txt --pheno-name $BN --1 \
    --covar pcs.eigenvec --covar-variance-standardize --glm hide-covar no-firth --out $O/m2b > /dev/null
  condense M2_b plink $O/m2b.*.glm.logistic
  $P2 --pfile base --extract tests.snplist --pheno pheno_b_iid.txt --pheno-name $BN --1 \
    --covar pcs.eigenvec --covar-variance-standardize --glm hide-covar firth-fallback --out $O/m2fb > /dev/null
  condense M2f_b plink $O/m2fb.*.glm.logistic.hybrid
  $RG --step 1 --bt --pgen rg --extract prune.prune.in --phenoFile pheno_b.txt --phenoColList $BN \
    --covarFile pcs_fidiid.txt --bsize 1000 --lowmem --lowmem-prefix $O/tmpb --out $O/s1b > /dev/null
  stamp "REGENIE step 1 binary"
  $RG --step 2 --bt --firth --approx --pgen rg --extract tests.snplist --phenoFile pheno_b.txt \
    --phenoColList $BN --covarFile pcs_fidiid.txt --pred $O/s1b_pred.list --bsize 400 --out $O/m3b > /dev/null
  condense M3_b regenie $O/m3b_*.regenie
  $RG --step 2 --bt --spa --pgen rg --extract tests.snplist --phenoFile pheno_b.txt \
    --phenoColList $BN --covarFile pcs_fidiid.txt --pred $O/s1b_pred.list --bsize 400 --out $O/m3sb > /dev/null
  condense M3s_b regenie $O/m3sb_*.regenie
}

case $PART in
  setup) setup ;;
  q) rm -f $RES; setup; run_q; cp keep200.raw $INPUTS/ ;;
  b) run_b ;;
  all) rm -f $RES; setup; run_q; run_b; cp keep200.raw $INPUTS/ ;;
esac
rm -rf $O
stamp "$PART done"
