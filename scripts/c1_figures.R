# figures: type I error, MR forest, PCSK9 locus
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
dir.create("figures", showWarnings = FALSE)
okabe <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#999999")
theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid.minor = element_blank(), legend.title = element_blank()))
save <- function(p, f, w = 6.5, h = 3.4) ggsave(file.path("figures", f), p, width = w, height = h, dpi = 150, bg = "white")

t1 <- fread("results/type1_error.csv")[method != "M1"]
lab <- c(M2 = "10 PCs", Mpop = "Population labels (post hoc)", M3 = "REGENIE + 10 PCs",
         M2f = "10 PCs, Firth fallback", M3s = "REGENIE + SPA")
t1[, method_lab := factor(lab[method], unique(lab))]
t1[, maf := factor(ifelse(maf_bin == "low", "MAF 0.01-0.05", "MAF >= 0.05"), c("MAF >= 0.05", "MAF 0.01-0.05"))]
t1[, kind := factor(kind, c("quantitative", "B20", "B05", "B02"),
                    c("Quantitative", "Binary, 20% cases", "Binary, 5% cases", "Binary, 2% cases"))]
band <- data.table(kind = factor(c("Quantitative", "Binary, 20% cases", "Binary, 5% cases", "Binary, 2% cases"),
                                 levels(t1$kind)), lo = c(0.8, 0.5, 0.5, 0.5), hi = c(1.25, 2, 2, 2))
band <- band[kind %in% t1$kind]
p1 <- ggplot(t1, aes(y = method_lab, colour = maf)) +
  geom_rect(data = band, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
            fill = "grey90") +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_pointrange(aes(x = `ratio_1e-4`, xmin = lo, xmax = hi), position = position_dodge(width = 0.5), size = 0.3) +
  scale_colour_manual(values = okabe[c(1, 6)]) +
  scale_x_log10(breaks = c(0.5, 0.8, 1, 1.25, 1.5, 2, 3)) +
  facet_wrap(~kind, ncol = 2) +
  labs(x = "Type I error at 1e-4 / nominal (95% bootstrap interval); grey = calibrated band", y = NULL,
       caption = "Null SNPs on chromosomes 12-22, 50 replicates per panel. No adjustment (M1) is off the scale (ratio above 1,000).") +
  theme(legend.position = "top")
save(p1, "type1_error.png", h = if (length(unique(t1$kind)) > 1) 5 else 2.8)

mr <- fread("results/mr_results.csv")[method == "IVW"]
mr[, outcome := factor(outcome, c("ihd", "t2d"), c("Ischemic heart disease", "Type 2 diabetes"))]
mr[, gene := factor(gene, rev(c("PCSK9", "HMGCR", "NPC1L1", "ANGPTL3")))]
p2 <- ggplot(mr, aes(OR_per_SD_lower, gene)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = okabe[1], size = 0.3) +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f), %d SNPs", OR_per_SD_lower, lo, hi, n_snps)),
            nudge_y = 0.32, size = 2.6, colour = "grey25") +
  scale_x_log10(breaks = c(0.3, 0.5, 1, 2), limits = c(0.28, 3)) + facet_wrap(~outcome) +
  labs(x = "Odds ratio per 1 SD genetically lower LDL (IVW, 95% CI)", y = NULL)
save(p2, "mr_forest.png", h = 3)

rd <- function(t) {
  d <- fread(sprintf("data/panukb/PCSK9_%s.tsv", t))[low_confidence_EUR == FALSE]
  d[, .(pos, lp = neglog10_pval_EUR, trait = t)]
}
loc <- rbind(rd("ldl"), rd("ihd"))
loc[, trait := factor(trait, c("ldl", "ihd"), c("LDL cholesterol", "Ischemic heart disease"))]
g <- fread("inputs/genes.csv")[gene == "PCSK9"]
p3 <- ggplot(loc, aes(pos / 1e6, lp)) +
  annotate("rect", xmin = g$start / 1e6, xmax = g$end / 1e6, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_point(size = 0.5, colour = okabe[1], alpha = 0.7) +
  facet_wrap(~trait, ncol = 1, scales = "free_y") +
  labs(x = "Chromosome 1 position (Mb, GRCh37); grey = PCSK9", y = "-log10 p (Pan-UKB, Europeans)")
save(p3, "pcsk9_locus.png", h = 3.6)
cat("figures written\n")
