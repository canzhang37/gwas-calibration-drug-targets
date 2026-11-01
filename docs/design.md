# Design: GWAS calibration under population structure, and genetic support for lipid drug targets

Written before any genotype was downloaded and before any summary statistic was read (only the Pan-UKB phenotype manifest was opened, to find file names and case counts). Changes after this point go into the amendments section at the end with the reason.

## Questions

1. In a sample with strong population structure, which association method keeps type I error at its nominal level for null variants: plain regression, regression on principal components (PCs), or a whole-genome regression mixed model (REGENIE)? Does that hold for binary traits with few cases and for low-frequency variants, where Firth and saddlepoint (SPA) corrections are meant to help? And at this sample size, does the mixed model gain power over PCs?
2. Using public UK Biobank summary statistics, do variants near the genes of approved LDL-lowering drugs (HMGCR, PCSK9, NPC1L1) and of ANGPTL3 predict lower ischemic heart disease risk, as the drugs do, and does the known on-target diabetes risk of HMGCR inhibition (statins) show up in the genetics?

## Part A: simulation on real genotypes

### Data

1000 Genomes phase 3, primary release (build 37, 2,504 people from 26 populations in 5 super-populations), in the PLINK 2 format distributed on the PLINK 2 resources page. The genotypes are real; only the traits are simulated, so the truth is known.

### Variant sets

- Autosomal biallelic SNPs with no missing calls.
- Step-1 / PCA set: minor allele frequency (MAF) at least 0.05, LD-pruned with PLINK 2 (`--indep-pairwise 1000kb 1 0.2`), all autosomes.
- Causal pool: MAF at least 0.05 on chromosomes 1-11, a random 20,000 SNPs (seed 2027).
- Test set: SNPs on chromosomes 12-22 with MAF at least 0.01, a random 300,000 (seed 2027), reported in two bins: MAF 0.01-0.05 and at least 0.05. No causal variant lies on chromosomes 12-22, so every test-set association is a false positive.
- PCs: 10 PCs from the pruned set (PLINK 2 `--pca`).
- Relatedness is not removed (1000 Genomes keeps a few related pairs); the number of pairs with KING kinship above 0.0884 is reported.

### Trait model (quantitative)

For each replicate, independently:

- Polygenic part: 500 causal SNPs drawn from the causal pool, effects on standardized genotypes drawn from N(0, 1) and scaled so this part has variance 0.25.
- Large-effect part: 10 more causal SNPs from the pool, 5 explaining 1% and 5 explaining 2% of trait variance (variance of the standardized genotype times the squared effect).
- Population environment: 26 population means drawn from N(0, 1), scaled to variance 0.10. This is ancestry-linked environment that genotypes do not cause.
- Noise: normal, variance chosen so the total variance is 1 (0.50 in expectation).

50 replicates. Binary traits threshold the same 50 liabilities at prevalence K = 0.20, 0.05 and 0.02 (about 500, 125 and 50 cases), giving 150 binary traits.

### Methods

| ID | Method | Quantitative | Binary |
|---|---|---|---|
| M1 | No adjustment | PLINK 2 `--glm` linear | PLINK 2 logistic, `no-firth` |
| M2 | 10 PCs | PLINK 2 linear + PCs | PLINK 2 logistic + PCs, `no-firth` |
| M2f | 10 PCs, Firth fallback | - | PLINK 2 logistic + PCs, `firth-fallback` |
| M3 | REGENIE + 10 PCs | step 1 ridge on the pruned set, LOCO; step 2 linear | step 2 with `--firth --approx` |
| M3s | REGENIE + 10 PCs, SPA | - | step 2 with `--spa` |

A test that returns no p-value (failed fit) counts as not rejected; the share of such tests is reported.

### Outcomes

- Type I error at alpha = 1e-4 over test-set SNPs, pooled over replicates, as a ratio to nominal, with a 95% interval from a bootstrap over replicates. Also reported at 1e-6 and 5e-8 (counts only; too few expected to estimate a ratio).
- Genomic inflation factor (lambda GC, median chi-square / 0.455) per trait, summarized by the median over replicates.
- Power at 5e-8 and mean chi-square at the 10 large-effect causal SNPs (quantitative traits).
- "Calibrated" means the 1e-4 ratio lies in [0.8, 1.25] for quantitative traits and in [0.5, 2] for binary traits (fewer tests reach small p-values when cases are few).

### Predictions

- P1. M1 on quantitative traits: type I error ratio above 2 in both MAF bins, and median lambda GC above 1.3.
- P2. M2 on quantitative traits: ratio above 1.25 in the MAF at least 0.05 bin (population-level environment that 10 PCs do not absorb).
- P3. M3 on quantitative traits: ratio in [0.8, 1.25] in both MAF bins.
- P4. Binary traits at K = 0.02, MAF 0.01-0.05: M2 ratio outside [0.5, 2]; M3 and M3s inside [0.5, 2].
- P5. Power: M3 mean chi-square at the large-effect SNPs is between 0.95 and 1.10 times M2's. At n = 2,504 a ridge prediction of a polygenic trait from about 10^5 SNPs carries little signal, so the mixed model should not gain much power here (its power gains come at biobank scale).

### Checks (must pass before results are reported)

- K1. The simulator's realized variance shares (polygenic, large-effect, environment, noise), averaged over replicates, are each within 0.02 of their targets.
- K2. M1 and M2 p-values from PLINK 2 match an independent least-squares fit in Python for 200 random SNP-trait pairs (largest absolute difference in log10 p below 0.01).
- K3. REGENIE step 2 run without step-1 predictions (`--ignore-pred`) agrees with M2 on quantitative traits (correlation of -log10 p above 0.999), so any M3-M2 difference comes from step 1.

## Part B: drug-target genetics from UK Biobank summary statistics

### Data

Pan-UK Biobank (Pan-UKB) summary statistics, European-ancestry columns, read by genomic range through the tabix index without downloading whole files:

- Exposure: LDL direct (biomarker 30780, inverse-rank-normal transformed, so effects are per SD), about 400,000 people.
- Primary outcome: ischemic heart disease (phecode 411; 37,672 cases, 382,052 controls).
- Secondary outcome: type 2 diabetes (phecode 250.2; 22,768 cases, 396,181 controls).
- LD reference: the 503 European-ancestry people in 1000 Genomes phase 3 (Part A data).
- Gene coordinates: Ensembl GRCh37 (HMGCR chr5:74,632,154-74,657,929; PCSK9 chr1:55,505,221-55,530,525; NPC1L1 chr7:44,552,134-44,580,914; ANGPTL3 chr1:63,063,158-63,071,830).

Exposure and outcomes come from the same cohort (sample overlap). With strong instruments (F statistic well above 10) the resulting bias toward the observational association is small; F statistics are reported.

### Instruments

Variants within 100 kb of the gene, LDL p < 5e-8, MAF at least 0.01 in 1000 Genomes Europeans, not flagged low confidence, alleles matched to 1000 Genomes (A/T and C/G variants with MAF above 0.42 dropped), clumped at r^2 < 0.1 in 1000 Genomes Europeans keeping the smallest LDL p. Effects are aligned to the LDL-lowering allele.

### Analysis

- Inverse-variance weighted (IVW) Mendelian randomization per gene, reported as the odds ratio per 1 SD lower LDL. Sensitivity: MR-Egger and weighted median when a gene has at least 3 instruments, and correlated-instrument IVW with the 1000 Genomes LD matrix (r^2 < 0.4).
- Colocalization (coloc, approximate Bayes factors, default priors) of LDL with each outcome in each window.

### Predictions

- P6. Ischemic heart disease: the IVW odds ratio per SD lower LDL is below 1 for HMGCR, PCSK9 and NPC1L1, with the 95% interval below 1 for at least two of the three.
- P7. ANGPTL3: point estimate below 1.
- P8. Type 2 diabetes: the HMGCR odds ratio per SD lower LDL is above 1 (the statin on-target effect).
- P9. Colocalization of LDL and ischemic heart disease: posterior probability of a shared variant (PP.H4) at least 0.8 at HMGCR and PCSK9.

### Checks

- K4. Allele frequencies of the instruments in Pan-UKB Europeans and 1000 Genomes Europeans correlate above 0.95 after harmonization.
- K5. The range reader returns rows sorted by position, all inside the requested range, and for one test window the same rows as a sequential scan from the first indexed block to the end of the window.

## Software

PLINK 2 and REGENIE (Linux builds, run in WSL Ubuntu), Python (NumPy, pandas, SciPy) for simulation and checks, R (MendelianRandomization, coloc, ggplot2) for Part B and figures.

## Amendments

1. Post hoc, after P3 failed on the quantitative traits. Under M2 and M3 the null SNPs stayed inflated (type I error ratio 1.28-1.68). To see whether the population-level environment in the trait model causes this, the quantitative traits were also tested with 25 population indicators as covariates instead of PCs (`scripts/a7_posthoc_population.sh`; PCs and indicators together are collinear in 1000 Genomes). Labeled post hoc wherever it appears.

Implementation notes (no change to the design): the causal pool is drawn after the MAF filter (PLINK 2 applies `--thin-count` before `--maf`, which first gave 11,409 SNPs instead of 20,000); the REGENIE MKL build is used because the default build ran single-threaded; the 150 binary traits run on MSI (`msi/`), since the laptop would need three to four hours, and the bundle checks that MSI rebuilds the same variant file and reproduces two laptop traits before the full run.
