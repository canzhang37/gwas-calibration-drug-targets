# Part B: drug-target MR and coloc, four LDL genes
suppressPackageStartupMessages({
  library(data.table); library(MendelianRandomization); library(coloc)
})
genes <- fread("inputs/genes.csv")
n_case <- list(ihd = c(37672, 382052), t2d = c(22768, 396181))
N_LDL <- 400223
P_GW <- -log10(5e-8)

read_pan <- function(g, t) {
  d <- fread(sprintf("data/panukb/%s_%s.tsv", g, t))
  d[, id := paste(chr, pos, ref, alt, sep = ":")]
  d[low_confidence_EUR == FALSE & !is.na(beta_EUR) & !is.na(se_EUR)]
}

read_ld <- function(g) {
  raw <- fread(sprintf("data/ld/%s.raw", g))
  G <- as.matrix(raw[, -(1:6)])
  colnames(G) <- sub("_[ACGT]+$", "", colnames(G))
  fq <- fread(sprintf("data/ld/%s.afreq", g))
  list(G = G, af = setNames(fq$ALT_FREQS, fq$ID))
}

clump <- function(d, G, r2) {
  d <- d[order(-neglog10_pval_EUR)]
  keep <- character(0)
  for (s in d$id) {
    if (!length(keep) || all(cor(G[, s], G[, keep, drop = FALSE])^2 < r2)) keep <- c(keep, s)
  }
  d[id %in% keep]
}

# single instrument: Wald ratio, first-order SE
wald <- function(bx, by, byse) {
  est <- by / bx; se <- byse / abs(bx)
  list(Estimate = est, CILower = est - 1.96 * se, CIUpper = est + 1.96 * se,
       Pvalue = 2 * pnorm(-abs(est / se)))
}

palin <- function(a, b) paste0(a, b) %in% c("AT", "TA", "CG", "GC")

mr_rows <- list(); inst_rows <- list(); coloc_rows <- list(); k4 <- list()
for (i in seq_len(nrow(genes))) {
  g <- genes$gene[i]
  ld <- read_ld(g)
  x <- read_pan(g, "ldl")[id %in% colnames(ld$G)]
  x[, af_1kg := ld$af[id]][, maf_1kg := pmin(af_1kg, 1 - af_1kg)]
  x <- x[maf_1kg >= 0.01 & !(palin(ref, alt) & maf_1kg > 0.42)]
  sig <- x[neglog10_pval_EUR > P_GW]
  k4[[g]] <- data.table(gene = g, n = nrow(sig), r = cor(sig$af_EUR, sig$af_1kg))

  for (out in c("ihd", "t2d")) {
    y <- read_pan(g, out)[, .(id, by = beta_EUR, byse = se_EUR, af_y = af_controls_EUR)]
    for (thr in c(0.1, 0.4)) {
      iv <- merge(clump(sig[id %in% y$id], ld$G, thr), y, by = "id")
      # align each instrument to its LDL-lowering allele
      s <- -sign(iv$beta_EUR)
      bx <- iv$beta_EUR * s; by <- iv$by * s
      rho <- cor(ld$G[, iv$id, drop = FALSE]) * outer(s, s)
      F <- mean((bx / iv$se_EUR)^2)
      if (nrow(iv) < 2 && thr == 0.4) next
      if (thr == 0.1) {
        inp <- mr_input(bx = bx, bxse = iv$se_EUR, by = by, byse = iv$byse, snps = iv$id)
        fits <- list(IVW = if (nrow(iv) > 1) mr_ivw(inp) else wald(bx, by, iv$byse))
        if (nrow(iv) >= 3) {
          fits$Egger <- mr_egger(inp)
          fits$`Weighted median` <- mr_median(inp, iterations = 10000, seed = 2027)
        }
        # fixed-effect IVW by hand, checks the package
        w <- 1 / iv$byse^2
        hand <- sum(bx * by * w) / sum(bx^2 * w)
        pkg <- if (nrow(iv) > 1) mr_ivw(inp, model = "fixed")$Estimate else hand
        inst_rows[[paste(g, out)]] <- data.table(gene = g, outcome = out, id = iv$id,
          bx = bx, bxse = iv$se_EUR, by = by, byse = iv$byse, maf_1kg = iv$maf_1kg)
      } else {
        inp <- mr_input(bx = bx, bxse = iv$se_EUR, by = by, byse = iv$byse,
                        correlation = rho, snps = iv$id)
        fits <- list(`IVW correlated (r2 < 0.4)` = mr_ivw(inp, correl = TRUE))
        hand <- pkg <- NA
      }
      for (m in names(fits)) {
        f <- fits[[m]]
        est <- if (m == "Egger") f$Estimate else f$Estimate
        lo <- if (m == "Egger") f$CILower.Est else f$CILower
        hi <- if (m == "Egger") f$CIUpper.Est else f$CIUpper
        p <- if (m == "Egger") f$Pvalue.Est else f$Pvalue
        # theta per SD higher LDL; reported per SD lower
        mr_rows[[length(mr_rows) + 1]] <- data.table(gene = g, outcome = out, method = m,
          n_snps = nrow(iv), F_mean = F, OR_per_SD_lower = exp(-est),
          lo = exp(-hi), hi = exp(-lo), p = p, hand_minus_pkg = hand - pkg)
      }
    }

    # colocalization over the whole window
    z <- merge(x[, .(id, bx = beta_EUR, vx = se_EUR^2, maf_1kg)],
               read_pan(g, out)[, .(id, b2 = beta_EUR, v2 = se_EUR^2)], by = "id")
    nc <- n_case[[out]]
    cc <- coloc.abf(
      dataset1 = list(beta = z$bx, varbeta = z$vx, snp = z$id, type = "quant", sdY = 1, N = N_LDL),
      dataset2 = list(beta = z$b2, varbeta = z$v2, snp = z$id, type = "cc",
                      s = nc[1] / sum(nc), N = sum(nc)))
    pp <- cc$summary
    coloc_rows[[paste(g, out)]] <- data.table(gene = g, outcome = out, nsnps = pp[["nsnps"]],
      PP.H0 = pp[["PP.H0.abf"]], PP.H1 = pp[["PP.H1.abf"]], PP.H2 = pp[["PP.H2.abf"]],
      PP.H3 = pp[["PP.H3.abf"]], PP.H4 = pp[["PP.H4.abf"]],
      lead_outcome_p = 10^-max(read_pan(g, out)$neglog10_pval_EUR))
  }
}
dir.create("results", showWarnings = FALSE)
mr <- rbindlist(mr_rows); fwrite(mr, "results/mr_results.csv")
fwrite(rbindlist(inst_rows), "results/mr_instruments.csv")
co <- rbindlist(coloc_rows); fwrite(co, "results/coloc.csv")
k4 <- rbindlist(k4); fwrite(k4, "results/k4_allele_freq.csv")
print(k4)
cat("K4", if (all(k4$r > 0.95)) "PASS" else "FAIL", "\n")
cat("IVW hand vs package, largest gap:", max(abs(mr$hand_minus_pkg), na.rm = TRUE), "\n")
print(mr[, .(gene, outcome, method, n_snps, F = round(F_mean), OR = round(OR_per_SD_lower, 3),
             lo = round(lo, 3), hi = round(hi, 3), p = signif(p, 2))])
print(co[, .(gene, outcome, nsnps, H3 = round(PP.H3, 3), H4 = round(PP.H4, 3))])
