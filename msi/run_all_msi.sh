#!/bin/bash
# Part A binary traits on MSI (150 x 5 methods)
# usage and settings: msi/README.md
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
MODE=${1:-run}
USER=${USER:-$(whoami)}
EMAIL=${EMAIL:-${USER}@umn.edu}
ACCOUNT=${ACCOUNT:-}
PARTITION=${PARTITION:-msismall}
DRYRUN=${DRYRUN:-0}
NJOBS=$(grep -c . msi/jobs.txt)
D=$ROOT/o_data
die() { echo "ERROR: $*" >&2; exit 1; }

if [ "$MODE" = "setup" ]; then
  if [ -n "${PY_MODULE:-}" ]; then module load "$PY_MODULE" || die "cannot load $PY_MODULE"
  else module load python3 2>/dev/null || true; fi
  python3 -c 'import sys; assert sys.version_info[:2] >= (3, 8), sys.version' \
    || die "need Python 3.8 or later; set PY_MODULE (see: module avail python)"
  mkdir -p $D/bin/mkl && cd $D/bin
  wget -q -O plink2.zip https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_avx2_20261024.zip
  wget -q -O regenie_mkl.zip https://github.com/rgcgithub/regenie/releases/download/v4.1.3/regenie_v4.1.3.gz_x86_64_Linux_mkl.zip
  python3 -m zipfile -e plink2.zip . && python3 -m zipfile -e regenie_mkl.zip mkl/
  chmod +x plink2 mkl/regenie_v4.1.3.gz_x86_64_Linux_mkl
  ./plink2 --version
  ldd mkl/regenie_v4.1.3.gz_x86_64_Linux_mkl | grep "not found" \
    && die "REGENIE is missing a library (often libgomp): try 'module load gcc' and rerun setup"
  mkl/regenie_v4.1.3.gz_x86_64_Linux_mkl --version
  echo "setup done"
  exit 0
fi

[ -x $D/bin/plink2 ] || die "run 'bash msi/run_all_msi.sh setup' first"

if [ "$MODE" = "data" ]; then
  cd $D
  wget -q -O all_phase3.psam "https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1"
  wget -q -O all_phase3.pvar.zst "https://www.dropbox.com/s/c95n8quqwqww4s0/all_phase3_noannot.pvar.zst?dl=1"
  wget -q -O all_phase3.pgen.zst "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1"
  ls -l $D
  echo "data done"
  exit 0
fi

if [ "$MODE" = "status" ]; then
  [ -f $D/work/base.pgen ] && echo "prep: base genotypes ready" || echo "prep: not done"
  [ -f output_smoke/smoke.txt ] && cat output_smoke/smoke.txt
  echo "finished tasks: $(ls output/task_*.done 2>/dev/null | wc -l)/$NJOBS"
  squeue -u "$USER" -o "%.12i %.18j %.8T %.10M %.20R" || true
  for f in $(ls -t logs/*.out 2>/dev/null | head -3); do echo "== $f"; tail -3 "$f"; done
  exit 0
fi

[ -f $D/all_phase3.pgen.zst ] || die "genotypes missing; run 'bash msi/run_all_msi.sh data' first"
mkdir -p logs
COMMON=(--chdir="$ROOT" --mail-user="$EMAIL" --mail-type=FAIL --kill-on-invalid-dep=yes --parsable --partition="$PARTITION")
if [ -n "$ACCOUNT" ]; then COMMON+=(--account="$ACCOUNT"); fi

submit() {  # prints the job id
  if [ "$DRYRUN" = "1" ]; then echo "sbatch $*" >&2; echo "DRY$RANDOM"; else sbatch "$@"; fi
}
step() {  # step name, then sbatch options
  local s=$1; shift
  submit "${COMMON[@]}" "$@" --job-name="gwas_$s" --export="ALL,ROOT=$ROOT,STEP=$s" msi/slurm_step.sh
}
array_job() {  # optional dependency
  step job ${1:+--dependency=afterok:"$1"} --cpus-per-task=8 --mem=16G --time=01:30:00 \
    --array=1-"$NJOBS" --output=logs/task_%A_%a.out
}
merge_job() {  # after all array tasks, failed included
  step merge --dependency=afterany:"$1" --cpus-per-task=1 --mem=2G --time=00:15:00 \
    --mail-type=END,FAIL --output=logs/merge_%j.out
}

case "$MODE" in
  run)
    P=$(step prep --cpus-per-task=8 --mem=16G --time=01:00:00 --output=logs/prep_%j.out)
    S=$(step smoke --dependency=afterok:"$P" --cpus-per-task=8 --mem=16G --time=00:45:00 --output=logs/smoke_%j.out)
    A=$(array_job "$S")
    M=$(merge_job "$A")
    echo "prep $P -> smoke $S -> array $A ($NJOBS tasks) -> merge $M"
    ;;
  rerun)
    A=$(array_job "")
    M=$(merge_job "$A")
    echo "array $A ($NJOBS tasks) -> merge $M"
    ;;
  merge)
    M=$(step merge --cpus-per-task=1 --mem=2G --time=00:15:00 --output=logs/merge_%j.out)
    echo "merge $M"
    ;;
  *)
    die "unknown mode $MODE"
    ;;
esac
echo "check with: bash msi/run_all_msi.sh status"
