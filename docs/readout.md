# Readout

Predictions and checks are those in `docs/design.md`, written before the data were opened. Post hoc results are labeled.

## Checks

| Check | Result |
|---|---|
| K1 simulator variance shares | Pass: mean shares 0.252 polygenic, 0.148 large-effect, 0.101 environment, 0.504 noise (largest gap 0.004) |
| K2 PLINK 2 vs least squares | Pass: largest log10 p gap 2e-6 (M1) and 1e-6 (M2) over 200 SNP-trait pairs |
| K3 REGENIE without step 1 vs PCs | Pass: correlation of -log10 p 1.000000 over 10,000 pairs |
| K4 allele frequencies, Pan-UKB vs 1000 Genomes Europeans | Pass: r = 0.991 to 0.998 by gene |
| K5 range reader | Pass: 2,499 rows in the PCSK9 window, sorted, in range, same as a sequential scan |

## Part A: calibration under population structure

1000 Genomes phase 3 after QC: 2,504 people, 12,123,596 biallelic SNPs with MAF at least 0.01; 283,759 SNPs after pruning for step 1 and the PCs; 42 pairs with kinship above 0.0884. The test set is 300,000 SNPs on chromosomes 12-22, where no causal variant lies.

Quantitative traits, 50 replicates, 15.0 million null tests:

| Method | Ratio at 1e-4, MAF >= 0.05 | Ratio at 1e-4, MAF 0.01-0.05 | Median lambda GC | Hits at 5e-8 (expected 0.75) |
|---|---|---|---|---|
| M1 no adjustment | 1,841 | 1,322 | 8.2 / 5.8 | 1,014,624 |
| M2 10 PCs | 1.50 (1.35-1.66) | 1.68 (1.50-1.88) | 1.04 / 1.06 | 5 |
| M3 REGENIE + 10 PCs | 1.28 (1.16-1.43) | 1.50 (1.31-1.71) | 1.02 / 1.04 | 1 |
| Population labels instead of PCs (post hoc) | 1.27 (1.16-1.38) | 1.44 (1.29-1.61) | 1.01 / 1.04 | 1 |

Intervals are 95% bootstrap intervals over replicates. Lambda GC is given for MAF at least 0.05 and MAF 0.01-0.05.

Power at the large-effect SNPs (250 per effect size):

| Method | Power at 5e-8, 1% variance | Power at 5e-8, 2% variance | Mean chi-square |
|---|---|---|---|
| M1 no adjustment | 0.42 | 0.74 | 45.9 |
| M2 10 PCs | 0.24 | 0.91 | 36.6 |
| M3 REGENIE + 10 PCs | 0.17 | 0.84 | 32.9 |
| Population labels (post hoc) | 0.30 | 0.94 | 39.1 |

M1's power is not usable: the same method rejects 1,841 times too often at null SNPs.

Binary traits (150 traits, P4): run on MSI, results pending.

| Prediction | Verdict |
|---|---|
| P1. M1 ratio above 2 and lambda above 1.3 | Met |
| P2. M2 ratio above 1.25 (MAF >= 0.05) | Met (1.50) |
| P3. M3 ratio in [0.8, 1.25], both MAF bins | Not met (1.28 and 1.50) |
| P4. Binary, 2% cases, low frequency: M2 outside [0.5, 2], M3 and M3s inside | Pending (MSI) |
| P5. M3 / M2 mean chi-square in [0.95, 1.10] | Not met (0.90) |

What this shows:

- Population structure across 26 populations makes unadjusted association tests useless: lambda GC 8.2, and about 1 in 11 null common SNPs reached p < 5e-8.
- Ten PCs remove most of it but leave 50% to 68% excess false positives at 1e-4. REGENIE removes about half of that excess for common variants (50% to 28%) and a quarter for low-frequency ones (68% to 50%). Neither is calibrated at this sample size, and low-frequency variants are worse than common ones under every method.
- Post hoc: replacing the PCs with the 26 population labels, which absorb the population-level environment exactly, leaves the same excess as REGENIE. So the remaining excess is not only that environment. A likely source, not tested here, is ancestry that varies within populations (the admixed American and African-American samples) and acts through the polygenic part of the trait.
- REGENIE did not gain power at n = 2,504; it lost some (mean chi-square 0.90 times that of PCs, prediction P5 expected 0.95 to 1.10). The step-1 predictor captures ancestry, and causal variants whose frequencies differ between populations lose part of their signal to it. Mixed-model power gains come at biobank sample sizes.

## Part B: genetic support for lipid drug targets

IVW Mendelian randomization with variants within 100 kb of each gene (LDL p < 5e-8, clumped at r^2 < 0.1), odds ratio per 1 SD genetically lower LDL:

| Gene (drug) | SNPs | Mean F | Ischemic heart disease | Type 2 diabetes |
|---|---|---|---|---|
| PCSK9 (alirocumab, evolocumab) | 18 | 212 | 0.52 (0.46-0.58) | 0.90 (0.76-1.07) |
| HMGCR (statins) | 11 | 144 | 0.82 (0.62-1.09) | 1.55 (1.21-1.98) |
| NPC1L1 (ezetimibe) | 2 | 137 | 0.45 (0.30-0.67) | 1.16 (0.70-1.93) |
| ANGPTL3 (evinacumab) | 3 | 129 | 1.36 (0.96-1.91) | 0.92 (0.60-1.41) |

Sensitivity analyses (`results/mr_results.csv`): for ischemic heart disease, the weighted median gives 0.55 for PCSK9 and 0.78 (0.64-0.96) for HMGCR, and correlated-instrument IVW with r^2 < 0.4 gives 0.56 for PCSK9, 0.81 for HMGCR and 0.43 for NPC1L1. The HMGCR diabetes estimate is 1.55 by weighted median and 1.60 by correlated IVW.

Colocalization of LDL with each outcome (PP.H4, probability of one shared causal variant):

| Gene | Ischemic heart disease | Type 2 diabetes |
|---|---|---|
| PCSK9 | 1.00 | 0.02 |
| HMGCR | 0.04 | 0.29 |
| NPC1L1 | 0.12 | 0.01 |
| ANGPTL3 | 0.02 | 0.01 |

| Prediction | Verdict |
|---|---|
| P6. IHD odds ratio below 1 for HMGCR, PCSK9 and NPC1L1, interval below 1 for at least two | Met (PCSK9 and NPC1L1 intervals below 1; HMGCR 0.82, interval crosses 1) |
| P7. ANGPTL3 point estimate below 1 | Not met (1.36) |
| P8. HMGCR type 2 diabetes odds ratio above 1 | Met (1.55) |
| P9. PP.H4 at least 0.8 at HMGCR and PCSK9 | Not met (PCSK9 1.00, HMGCR 0.04) |

What this shows:

- PCSK9 is the clearest case: 18 instruments, lower LDL predicts about half the odds of ischemic heart disease per SD, the LDL and heart disease signals share one variant (PP.H4 1.00), and no diabetes signal was detected. The top LDL variant in the window, chr1:55,505,647 G>T (the known PCSK9 R46L variant, 1.8% in Europeans; 0.42 SD lower LDL, p about 1e-448), is also the top heart disease variant (p about 3e-14).
- HMGCR reproduces the known statin trade-off: lower LDL through HMGCR predicts less heart disease (weaker here than for PCSK9; the IVW interval crosses 1, the weighted median does not) and more type 2 diabetes (1.55 per SD).
- The ANGPTL3 window does not support the target with common variants. The LDL-associated variants there do not share a signal with heart disease (PP.H4 0.02), and the point estimate goes the wrong way. ANGPTL3's protective evidence comes from rare loss-of-function carriers, and its main effect is on triglycerides; common LDL variants within 100 kb may act through neighboring genes. This is an interpretation after seeing the result.
- Exposure and outcomes come from the same UK Biobank participants. With mean F statistics of 129 to 212 the bias toward the observational association should be small, but it is not zero.

## Limitations

- 2,504 people is small for GWAS; the calibration results say what happens at this size and this degree of structure, not in a biobank.
- The simulated traits have one form of environment (population means) and a normal model; real traits can have other confounders.
- Part B uses summary statistics from one cohort with European ancestry only, and the 503 1000 Genomes Europeans as LD reference, which makes clumping and the LD matrix noisy.
- Gene windows of 100 kb can include variants that act through neighboring genes.
