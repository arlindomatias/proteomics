################################################################
# Hippocampal proteome after ACT chemotherapy (doxorubicin + cyclophosphamide + paclitaxel)
# McElroy et al. (2020) Int. J. Mol. Sci. -- supplementary file ijms-730741 (Table S1, raw iBAQ)
# Contrast: ACT vs Sham (n = 4 per group). Positive logFC = higher in ACT.
# Input: processed.csv (built by prepare.sh). Run from this folder.
################################################################

######################## Packages ########################
suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(DOSE)
})

######################## Adjusts ########################
set.seed(42)                      # imputation is stochastic
options(scipen = 0)
palette(c("#3399FF", "#CC0000", "#66A61E"))
group_colors <- c(Sham = "#3399FF", ACT = "#CC0000")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Save a table as CSV, adding the extension if missing
save <- function(table, file_name) {
  if (!grepl("\\.csv$", file_name, ignore.case = TRUE)) file_name <- paste0(file_name, ".csv")
  write.csv(table, file = file.path("results", file_name), row.names = FALSE)
}

# Save a ggplot / base plot to figures/
save_plot <- function(plot, file_name, width = 7, height = 6) {
  ggsave(file.path("figures", file_name), plot, width = width, height = height)
}

# Thresholds for calling a protein differentially abundant:
# BH FDR < 0.05 first; if nothing passes, fall back to unadjusted P < 0.01.
fdr_cutoff     <- 0.05
nominal_cutoff <- 0.01
fc_cutoff      <- 0.5    # |log2FC|

######################## Import Data ########################
raw <- read.csv("processed.csv", check.names = FALSE)

sample_cols <- grep("^(Sham|ACT)_", colnames(raw), value = TRUE)
annotation <- raw %>%
  transmute(Protein,
            protein_name = `Protein names`,
            SYMBOL = sub(";.*", "", `Gene names`))   # first symbol of the protein group

intensity <- as.matrix(raw[, sample_cols])
rownames(intensity) <- raw$Protein

# Sample metadata. Replicates 1-2 vs 3-4 differ strongly in depth (missing values),
# so the replicate pair is kept as a candidate batch factor (checked in the PCA below).
metadata <- data.frame(
  sample = sample_cols,
  group  = factor(sub("_.*", "", sample_cols), levels = c("Sham", "ACT")),
  batch  = factor(ifelse(as.integer(sub(".*_", "", sample_cols)) <= 2, "B1", "B2")),
  row.names = sample_cols
)

######################## Pre-processing ########################
## Missing values: iBAQ 0 = not quantified
intensity[intensity == 0] <- NA
log_int <- log2(intensity)

missing_per_sample <- colSums(is.na(log_int))
print(missing_per_sample)

## Filter: valid values in >= 3 of 4 replicates in at least one group
valid_by_group <- sapply(levels(metadata$group), function(g)
  rowSums(!is.na(log_int[, metadata$group == g, drop = FALSE])))
keep <- apply(valid_by_group >= 3, 1, any)
log_filt <- log_int[keep, ]
cat("Proteins kept after valid-value filter:", nrow(log_filt), "of", nrow(log_int), "\n")

## Normalization: median centering (shift every sample to the global median)
sample_medians <- apply(log_filt, 2, median, na.rm = TRUE)
log_norm <- sweep(log_filt, 2, sample_medians - median(sample_medians))

## Imputation: MinProb-style left-censored draws (missing values assumed MNAR, i.e. below
## detection). Per sample: N(mean - 1.8 * sd, 0.3 * sd), as in Perseus.
impute_minprob <- function(m, shift = 1.8, width = 0.3) {
  for (j in seq_len(ncol(m))) {
    miss <- is.na(m[, j])
    mu <- mean(m[, j], na.rm = TRUE); s <- sd(m[, j], na.rm = TRUE)
    m[miss, j] <- rnorm(sum(miss), mean = mu - shift * s, sd = width * s)
  }
  m
}
log_imp <- impute_minprob(log_norm)

## Densities before / after normalization and imputation
long <- function(m, stage) data.frame(m, check.names = FALSE) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "log2_iBAQ") %>%
  mutate(stage = stage)
dens <- bind_rows(long(log_filt, "1 filtered"), long(log_norm, "2 normalized"),
                  long(log_imp, "3 imputed")) %>%
  filter(!is.na(log2_iBAQ)) %>%
  mutate(group = sub("_.*", "", sample))
p_dens <- ggplot(dens, aes(log2_iBAQ, colour = sample, linetype = group)) +
  geom_density() + facet_wrap(~stage, ncol = 1) + theme_bw()
save_plot(p_dens, "density_preprocessing.pdf", height = 8)

## Boxplot after normalization (before imputation)
p_box <- long(log_norm, "normalized") %>% filter(!is.na(log2_iBAQ)) %>%
  mutate(group = sub("_.*", "", sample)) %>%
  ggplot(aes(sample, log2_iBAQ, fill = group)) + geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = group_colors) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_box, "boxplot_normalized.pdf")

## PCA of samples (imputed matrix)
pca <- prcomp(t(log_imp), scale. = FALSE)
var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- cbind(metadata, pca$x[, 1:3])
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = group, shape = batch, label = sample)) +
  geom_point(size = 4) + geom_text_repel(show.legend = FALSE) +
  scale_colour_manual(values = group_colors) +
  labs(x = paste0("PC1 (", var_expl[1], "%)"), y = paste0("PC2 (", var_expl[2], "%)")) +
  theme_bw()
save_plot(p_pca, "pca.pdf")
print(var_expl[1:4])

######################## Differential expression ########################
# Batch (replicate pair) is included as a blocking factor; see PCA for justification.
design <- model.matrix(~ 0 + group + batch, data = metadata)
colnames(design) <- sub("^group", "", colnames(design))
contrast <- makeContrasts(ACT_vs_Sham = ACT - Sham, levels = design)

fit <- lmFit(log_imp, design)
fit2 <- eBayes(contrasts.fit(fit, contrast), trend = TRUE)
tt <- topTable(fit2, coef = "ACT_vs_Sham", number = Inf, sort.by = "P")
tt$Protein <- rownames(tt)
tt <- left_join(tt, annotation, by = "Protein") %>%
  mutate(n_missing_Sham = valid_by_group[Protein, "Sham"] %>% {4 - .},
         n_missing_ACT  = valid_by_group[Protein, "ACT"]  %>% {4 - .}) %>%
  relocate(Protein, SYMBOL, protein_name)

# Significance column: FDR if anything passes it, otherwise unadjusted P < 0.01
if (any(tt$adj.P.Val < fdr_cutoff & abs(tt$logFC) > fc_cutoff)) {
  sig_col <- "adj.P.Val"; sig_cutoff <- fdr_cutoff
} else {
  sig_col <- "P.Value";   sig_cutoff <- nominal_cutoff
}
cat("Significance criterion:", sig_col, "<", sig_cutoff, "& |logFC| >", fc_cutoff, "\n")

tt$status <- with(tt, ifelse(tt[[sig_col]] < sig_cutoff & logFC >  fc_cutoff, "Up",
                      ifelse(tt[[sig_col]] < sig_cutoff & logFC < -fc_cutoff, "Down", "NS")))
tt_up   <- subset(tt, status == "Up")
tt_down <- subset(tt, status == "Down")
cat("Up in ACT:", nrow(tt_up), " Down in ACT:", nrow(tt_down), "\n")

## Sensitivity analysis: same model on the normalized matrix WITHOUT imputation.
# limma drops NAs per protein; proteins absent from a whole group get no estimate (NA).
fit_noimp <- lmFit(log_norm, design)
fit2_noimp <- eBayes(contrasts.fit(fit_noimp, contrast), trend = TRUE)
tt_noimp <- topTable(fit2_noimp, coef = "ACT_vs_Sham", number = Inf, sort.by = "none")
tt_noimp$Protein <- rownames(tt_noimp)

tt <- tt %>%
  left_join(tt_noimp %>% dplyr::select(Protein, noimp_logFC = logFC, noimp_P.Value = P.Value,
                                noimp_adj.P.Val = adj.P.Val), by = "Protein") %>%
  mutate(noimp_sig = if (sig_col == "adj.P.Val") noimp_adj.P.Val else noimp_P.Value,
         # DE call holds without imputation: same direction, same cutoffs
         robust_noimp = status != "NS" & !is.na(noimp_sig) & noimp_sig < sig_cutoff &
                        abs(noimp_logFC) > fc_cutoff & sign(noimp_logFC) == sign(logFC)) %>%
  dplyr::select(-noimp_sig)
tt_up   <- subset(tt, status == "Up")
tt_down <- subset(tt, status == "Down")

cat("Correlation of logFC, imputed vs not imputed:",
    round(cor(tt$logFC, tt$noimp_logFC, use = "complete.obs"), 3), "\n")
cat("DE proteins robust without imputation:", sum(tt$robust_noimp), "of", sum(tt$status != "NS"), "\n")
cat("No-imputation fit alone -", sig_col, "<", sig_cutoff, ":",
    sum(tt_noimp[[sig_col]] < sig_cutoff & abs(tt_noimp$logFC) > fc_cutoff, na.rm = TRUE), "proteins\n")

p_sens <- ggplot(tt, aes(logFC, noimp_logFC, colour = robust_noimp)) +
  geom_abline(colour = "grey") + geom_point(alpha = 0.4, size = 1) +
  geom_text_repel(data = subset(tt, status != "NS"), aes(label = SYMBOL), size = 2.5,
                  max.overlaps = 30, show.legend = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "grey50", `TRUE` = "#CC0000")) +
  labs(x = "log2FC (imputed)", y = "log2FC (no imputation)", colour = "DE robust\nwithout imputation") +
  theme_bw()
save_plot(p_sens, "sensitivity_imputation_logFC.pdf", width = 8, height = 7)

save(tt, "ACT_vs_Sham_limma_all")
save(rbind(tt_up, tt_down), "ACT_vs_Sham_limma_DE")

######################## Visualization ########################
## Volcano plot
p_volc <- EnhancedVolcano(tt,
  lab = tt$SYMBOL, x = "logFC", y = sig_col,
  pCutoff = sig_cutoff, FCcutoff = fc_cutoff,
  ylab = bquote(~-Log[10] ~ .(sig_col)),
  xlab = bquote(~Log[2] ~ "fold change (ACT / Sham)"),
  title = "Hippocampal proteome: ACT vs Sham",
  subtitle = paste0("limma; ", sig_col, " < ", sig_cutoff, ", |log2FC| > ", fc_cutoff),
  caption = paste("Total proteins:", nrow(tt)),
  pointSize = 2, labSize = 3, colAlpha = 0.5,
  col = c("gray", "#3399FF", "#009933", "#CC0000"),
  legendLabels = c("NS", "Log2FC", sig_col, paste(sig_col, "& Log2FC")),
  drawConnectors = TRUE, widthConnectors = 0.3, max.overlaps = 30)
save_plot(p_volc, "volcano.pdf", width = 9, height = 8)

## Heatmap of the differentially abundant proteins (row-scaled)
top_ids <- c(tt_up$Protein, tt_down$Protein)
hm <- log_imp[top_ids, ]
rownames(hm) <- make.unique(tt$SYMBOL[match(top_ids, tt$Protein)])
pdf("figures/heatmap_DE.pdf", width = 6, height = max(5, 0.15 * length(top_ids) + 2))
pheatmap(hm, scale = "row", cluster_cols = FALSE,
         annotation_col = metadata[, c("group", "batch")],
         annotation_colors = list(group = group_colors),
         main = paste0("DE proteins (", sig_col, " < ", sig_cutoff, ")"))
dev.off()

## P value / FDR histogram
pdf("figures/pvalue_histogram.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
hist(tt$P.Value, breaks = 50, col = "grey", main = "P value", xlab = "P")
hist(tt$adj.P.Val, breaks = 50, col = "grey", main = "FDR (BH)", xlab = "adj. P")
dev.off()

## QQ plot of moderated t
pdf("figures/qq_t.pdf", width = 5, height = 5)
qqt(fit2$t[, 1], df = fit2$df.total, main = "Moderated t vs theoretical")
abline(0, 1, col = "red")
dev.off()

## MD plot (logFC vs mean log intensity)
pdf("figures/md_plot.pdf", width = 6, height = 5)
plotMD(fit2, coef = 1, status = tt$status[match(rownames(fit2), tt$Protein)],
       main = "ACT vs Sham", values = c("Up", "Down"), col = c("#CC0000", "#3399FF"))
abline(h = 0, col = "grey")
dev.off()

## Mean-variance trend
pdf("figures/sa_plot.pdf", width = 6, height = 5)
plotSA(fit2, main = "Residual SD vs mean (trend = TRUE)")
dev.off()

######################## Enrichment ########################
# Universe = identified proteome after filtering, not the whole genome.
keytable <- bitr(unique(tt$SYMBOL[tt$SYMBOL != ""]), fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Mm.eg.db)
tt_entrez <- inner_join(tt, keytable, by = "SYMBOL") %>%
  group_by(ENTREZID) %>% slice_min(P.Value, n = 1, with_ties = FALSE) %>% ungroup()
universe <- tt_entrez$ENTREZID

## Ranked list for GSEA: signed -log10(P) keeps both direction and evidence
gene_list <- with(tt_entrez, sign(logFC) * -log10(P.Value))
names(gene_list) <- tt_entrez$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

## GSEA on GO (permutation-based: reseed so results do not depend on earlier random draws)
set.seed(42)
gse_go <- gseGO(geneList = gene_list, ont = "ALL", keyType = "ENTREZID",
                OrgDb = org.Mm.eg.db, minGSSize = 10, maxGSSize = 500,
                pvalueCutoff = 0.05, pAdjustMethod = "BH", eps = 0, verbose = FALSE)
gse_go <- setReadable(gse_go, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
save(as.data.frame(gse_go), "gsea_GO")
cat("GSEA GO terms (FDR < 0.05):", nrow(gse_go), "\n")

if (nrow(gse_go) > 0) {
  save_plot(dotplot(gse_go, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign),
            "gsea_GO_dotplot.pdf", width = 11, height = 8)
  save_plot(ridgeplot(gse_go, showCategory = 15) + labs(x = "signed -log10(P)"),
            "gsea_GO_ridgeplot.pdf", width = 9, height = 8)
  gse_go_sim <- pairwise_termsim(gse_go)
  save_plot(emapplot(gse_go_sim, showCategory = 30), "gsea_GO_emapplot.pdf", width = 10, height = 9)
  save_plot(gseaplot2(gse_go, geneSetID = 1, title = gse_go$Description[1], pvalue_table = TRUE),
            "gsea_GO_top1.pdf", width = 8, height = 6)
}

## GSEA on KEGG (requires internet access to KEGG)
gse_kegg <- tryCatch(
  gseKEGG(geneList = gene_list, organism = "mmu", keyType = "ncbi-geneid",
          minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0, verbose = FALSE),
  error = function(e) { message("KEGG skipped: ", conditionMessage(e)); NULL })
if (!is.null(gse_kegg)) {
  gse_kegg <- setReadable(gse_kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  save(as.data.frame(gse_kegg), "gsea_KEGG")
  cat("GSEA KEGG pathways (FDR < 0.05):", nrow(gse_kegg), "\n")
  if (nrow(gse_kegg) > 0)
    save_plot(dotplot(gse_kegg, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign),
              "gsea_KEGG_dotplot.pdf", width = 11, height = 8)
}

## ORA of up / down proteins (same significance criterion as the DE call)
ora_sets <- list(
  up   = tt_entrez$ENTREZID[tt_entrez$status == "Up"],
  down = tt_entrez$ENTREZID[tt_entrez$status == "Down"])
ora <- compareCluster(ora_sets, fun = "enrichGO", OrgDb = org.Mm.eg.db, ont = "BP",
                      universe = universe, pvalueCutoff = 0.05, readable = TRUE)
if (!is.null(ora) && nrow(as.data.frame(ora)) > 0) {
  save(as.data.frame(ora), "ora_GO_BP_up_down")
  save_plot(dotplot(ora, showCategory = 10), "ora_GO_BP_dotplot.pdf", width = 9, height = 8)
}

sessionInfo()
