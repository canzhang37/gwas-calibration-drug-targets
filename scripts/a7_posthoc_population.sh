#!/usr/bin/env bash
# post hoc, amendment 1: 25 population labels
# residual inflation = population environment?
set -euo pipefail
PROJ=${PROJ:-$(cd "$(dirname "$0")/.." && pwd)}
D=$HOME/o_data
W=$D/work
P2="$D/bin/plink2 --memory 6500 --threads 16"
export O_WORK=$W
RES=$PROJ/data/work/posthoc.jsonl
rm -f $RES
cd $W
# indicators from the psam, first population = ref
# no PCs (collinear with the labels in 1000G)
python3 - <<'PY'
rows = [l.rstrip("\n").split("\t") for l in open("base.psam")]
ip = rows[0].index("Population")
pops = sorted({r[ip] for r in rows[1:]})[1:]
with open("pop_covar.txt", "w") as f:
    f.write("\t".join(["#IID"] + ["pop_" + p for p in pops]) + "\n")
    for r in rows[1:]:
        f.write("\t".join([r[0]] + ["1" if r[ip] == p else "0" for p in pops]) + "\n")
PY
QN=$(head -1 pheno_q.txt | cut -f3- | tr '\t' ',')
mkdir -p out_post
$P2 --pfile base --extract tests.snplist --pheno pheno_q_iid.txt --pheno-name $QN \
  --covar pop_covar.txt --glm hide-covar --out out_post/mpop > /dev/null
python3 $PROJ/scripts/a4_summarize.py Mpop_q plink $RES out_post/mpop.*.glm.linear > /dev/null
rm -rf out_post
echo "post hoc done: $(wc -l < $RES) traits"
