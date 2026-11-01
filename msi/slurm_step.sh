#!/bin/bash
# one step: prep, smoke, job or merge
set -euo pipefail
cd "${ROOT:?}"
if [ -n "${PY_MODULE:-}" ]; then module load "$PY_MODULE"; else module load python3 2>/dev/null || true; fi
export PROJ=$ROOT O_DATA=$ROOT/o_data INPUTS=$ROOT/msi/inputs
export PLINK2_BIN=$O_DATA/bin/plink2
export REGENIE_BIN=$O_DATA/bin/mkl/regenie_v4.1.3.gz_x86_64_Linux_mkl
export THREADS=${SLURM_CPUS_PER_TASK:-8} MEM=12000
W=$O_DATA/work

case "$STEP" in
  prep)
    BASE_ONLY=1 bash scripts/a1_prepare.sh
    # the laptop's variant set must be rebuilt exactly
    (cd $W && md5sum -c $INPUTS/base_pvar.md5)
    cp $INPUTS/pcs.eigenvec $INPUTS/prune.prune.in $INPUTS/testset.snplist $W/
    PART=setup bash scripts/a3_run_methods.sh
    echo "- end"
    ;;
  smoke)
    rm -rf output_smoke && mkdir -p output_smoke
    TRAITS=B20_01,B02_01 TAG=smoke RES=$ROOT/output_smoke/smoke.jsonl PART=b bash scripts/a3_run_methods.sh
    python3 msi/o_tools.py smoke | tee output_smoke/smoke.txt
    echo "- end"
    ;;
  job)
    TR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" msi/jobs.txt)
    echo "task $SLURM_ARRAY_TASK_ID: $TR"
    mkdir -p output
    # rerun starts the task afresh
    rm -f output/task_$SLURM_ARRAY_TASK_ID.jsonl output/task_$SLURM_ARRAY_TASK_ID.done
    TRAITS=$TR TAG=task$SLURM_ARRAY_TASK_ID RES=$ROOT/output/task_$SLURM_ARRAY_TASK_ID.jsonl PART=b \
      bash scripts/a3_run_methods.sh
    touch output/task_$SLURM_ARRAY_TASK_ID.done
    echo "- end"
    ;;
  merge)
    python3 msi/o_tools.py merge
    ;;
  *)
    echo "unknown STEP $STEP"; exit 1
    ;;
esac
