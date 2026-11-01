# Binary traits on MSI

Part A tests 150 simulated binary traits (prevalence 20%, 5% and 2%, 50 replicates each) on 300,000 null SNPs with five methods: PLINK 2 logistic regression without and with 10 PCs, PCs with Firth fallback, and REGENIE with Firth or saddlepoint corrections. On the laptop this takes three to four hours, so it runs here. The quantitative traits and Part B were run on the laptop.

## Steps

From a login node, in the unpacked bundle folder:

```
bash msi/run_all_msi.sh setup
bash msi/run_all_msi.sh data
bash msi/run_all_msi.sh
bash msi/run_all_msi.sh status
```

`setup` downloads PLINK 2 (alpha 7) and REGENIE 4.1.3 (MKL build) into `o_data/bin` and checks that REGENIE finds its libraries. If it reports a missing library, run `module load gcc` and `setup` again. It needs Python 3.8 or later (`module load python3`, or set `PY_MODULE`).

`data` downloads the 1000 Genomes phase 3 genotypes in PLINK 2 format (about 3 GB) from the PLINK 2 resources page into `o_data`.

The last command submits four steps, each started only if the one before it succeeded:

1. prep: rebuilds the quality-controlled genotypes (12.1 million SNPs) and checks that the variant file has the same MD5 sum as on the laptop. The pruned SNP set, PCs and test SNPs come from the laptop in `msi/inputs`.
2. smoke: runs all five methods on two traits and checks the counts against the laptop run in `msi/inputs/smoke_laptop.jsonl`.
3. an array of 15 tasks, 10 traits each (`msi/jobs.txt`).
4. merge.

## Resources

Each array task uses 8 cores and up to 16 GB for about 30 to 60 minutes, about 100 core-hours in all; prep and smoke take about 30 minutes each. `EMAIL`, `ACCOUNT`, `PARTITION` and `PY_MODULE` can be set before a command (see the top of `msi/run_all_msi.sh`).

## Results

The merge job writes `gwas_binary_results_<date>.zip` in the bundle root, with `output/binary_summary.jsonl` (one line per method and trait), `STATUS.txt` and the smoke check. Only this zip comes back to the laptop. If `STATUS.txt` lists unfinished tasks, `bash msi/run_all_msi.sh rerun` runs the array and the merge again.
