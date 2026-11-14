################################################################
# Hippocampal proteome after ACT chemotherapy (doxorubicin + cyclophosphamide + paclitaxel)
# McElroy et al. (2020) Int. J. Mol. Sci. -- supplementary file ijms-730741 (Table S1, raw iBAQ)
# Contrast: ACT vs Sham (n = 4 per group). Positive logFC = higher in ACT.
# Input: data/processed.csv (built by scripts/prepare.sh). Parameters: config/config.yaml.
# Usage: `Rscript scripts/<this>.R [config/config.yaml]` (or via `snakemake`). Paths are relative
# to the project root (proteomics/); for interactive use, set the working directory there.
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
  library(yaml)
})

######################## Adjusts ########################
# Move to the project root when run with Rscript (no-op in interactive sessions)
script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script_arg) == 1) setwd(file.path(dirname(normalizePath(script_arg)), ".."))

# Parameters: first command-line argument, default config/config.yaml
args <- commandArgs(trailingOnly = TRUE)
cfg <- read_yaml(if (length(args) >= 1) args[1] else "config/config.yaml")
results_dir <- cfg$results_dir
plots_dir   <- cfg$plots_dir

set.seed(cfg$seed)                # imputation is stochastic
options(scipen = 0)
palette(c("#3399FF", "#CC0000", "#66A61E"))
group_colors <- unlist(cfg$groups$colors)    # Sham blue, ACT red (as in bulk-rna-seq)

dir.create(results_dir, showWarnings = FALSE)   # tables
dir.create(plots_dir, showWarnings = FALSE)     # figures

# Save a table as CSV, adding the extension if missing
save <- function(table, file_name) {
  if (!grepl("\\.csv$", file_name, ignore.case = TRUE)) file_name <- paste0(file_name, ".csv")
  write.csv(table, file = file.path(results_dir, file_name), row.names = FALSE)
}

# Save a ggplot to the plots folder as PDF (for print) and PNG (embedded in the HTML report)
save_plot <- function(plot, file_name, width = 7, height = 6) {
  ggsave(file.path(plots_dir, file_name), plot, width = width, height = height)
  ggsave(file.path(plots_dir, sub("\\.pdf$", ".png", file_name)), plot,
         width = width, height = height, dpi = 150)
}

# Same for base-graphics plots: `code` is evaluated once per device
save_base <- function(file_name, width, height, code) {
  code <- substitute(code); env <- parent.frame()
  pdf(file.path(plots_dir, file_name), width = width, height = height)
  eval(code, env); invisible(dev.off())
  png(file.path(plots_dir, sub("\\.pdf$", ".png", file_name)),
      width = width, height = height, units = "in", res = 150)
  eval(code, env); invisible(dev.off())
}

# Thresholds for calling a protein differentially abundant:
# BH FDR first; if nothing passes, fall back to unadjusted P.
fdr_cutoff     <- cfg$de$fdr_cutoff
nominal_cutoff <- cfg$de$nominal_cutoff
fc_cutoff      <- cfg$de$logfc_cutoff    # |log2FC|

######################## Import Data ########################
raw <- read.csv(cfg$data$processed, check.names = FALSE)

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

## Filter: valid values in >= min_valid replicates in at least one group
valid_by_group <- sapply(levels(metadata$group), function(g)
  rowSums(!is.na(log_int[, metadata$group == g, drop = FALSE])))
keep <- apply(valid_by_group >= cfg$preprocessing$min_valid_per_group, 1, any)
log_filt <- log_int[keep, ]
cat("Proteins kept after valid-value filter:", nrow(log_filt), "of", nrow(log_int), "\n")

## Normalization: median centering (shift every sample to the global median)
sample_medians <- apply(log_filt, 2, median, na.rm = TRUE)
log_norm <- sweep(log_filt, 2, sample_medians - median(sample_medians))

## Imputation: MinProb-style left-censored draws (missing values assumed MNAR, i.e. below
## detection). Per sample: N(mean - shift * sd, width * sd); Perseus default 1.8 / 0.3.
impute_minprob <- function(m, shift = 1.8, width = 0.3) {
  for (j in seq_len(ncol(m))) {
    miss <- is.na(m[, j])
    mu <- mean(m[, j], na.rm = TRUE); s <- sd(m[, j], na.rm = TRUE)
    m[miss, j] <- rnorm(sum(miss), mean = mu - shift * s, sd = width * s)
  }
  m
}
log_imp <- impute_minprob(log_norm, shift = cfg$preprocessing$impute_shift,
                          width = cfg$preprocessing$impute_width)

# <figure name="density_preprocessing">
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
# </figure>

# <figure name="boxplot_normalized">
## Boxplot after normalization (before imputation)
p_box <- long(log_norm, "normalized") %>% filter(!is.na(log2_iBAQ)) %>%
  mutate(group = sub("_.*", "", sample)) %>%
  ggplot(aes(sample, log2_iBAQ, fill = group)) + geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = group_colors) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_box, "boxplot_normalized.pdf")
# </figure>

## Batch-corrected matrix, for visualization only (PCA, heatmap). The DE model does not use it:
## it keeps the raw matrix and estimates the batch as a term (`~ 0 + group + batch`).
## `design` protects the group effect from being removed together with the batch.
log_bc <- removeBatchEffect(log_imp, batch = metadata$batch,
                            design = model.matrix(~ group, data = metadata))

# <figure name="pca">
## PCA of samples (imputed, batch-corrected matrix)
pca <- prcomp(t(log_bc), scale. = FALSE)
var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- cbind(metadata, pca$x[, 1:3])
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = group, shape = batch, label = sample)) +
  geom_point(size = 4) + geom_text_repel(show.legend = FALSE) +
  scale_colour_manual(values = group_colors) +
  labs(x = paste0("PC1 (", var_expl[1], "%)"), y = paste0("PC2 (", var_expl[2], "%)"),
       title = "PCA, batch-corrected (removeBatchEffect)") +
  theme_bw()
save_plot(p_pca, "pca.pdf")
# </figure>
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
  mutate(n_missing_Sham = sum(metadata$group == "Sham") - valid_by_group[Protein, "Sham"],
         n_missing_ACT  = sum(metadata$group == "ACT")  - valid_by_group[Protein, "ACT"]) %>%
  relocate(Protein, SYMBOL, protein_name)

# Significance column: FDR if anything passes it, otherwise unadjusted P
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

# <figure name="sensitivity_imputation_logFC">
p_sens <- ggplot(tt, aes(logFC, noimp_logFC, colour = robust_noimp)) +
  geom_abline(colour = "grey") + geom_point(alpha = 0.4, size = 1) +
  geom_text_repel(data = subset(tt, status != "NS"), aes(label = SYMBOL), size = 2.5,
                  max.overlaps = 30, show.legend = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "grey50", `TRUE` = "#CC0000")) +
  labs(x = "log2FC (imputed)", y = "log2FC (no imputation)", colour = "DE robust\nwithout imputation") +
  theme_bw()
save_plot(p_sens, "sensitivity_imputation_logFC.pdf", width = 8, height = 7)
# </figure>

save(tt, "ACT_vs_Sham_limma_all")
save(rbind(tt_up, tt_down), "ACT_vs_Sham_limma_DE")

######################## Visualization ########################
# <figure name="volcano">
# Volcano plot: Up = red, Down = blue, NS = gray (same scheme as bulk-rna-seq)
keyvals <- ifelse(tt$status == "Up", "#CC0000", ifelse(tt$status == "Down", "#3399FF", "gray"))
names(keyvals) <- ifelse(keyvals == "#CC0000", "Up", ifelse(keyvals == "#3399FF", "Down", "NS"))
sig_label <- if (sig_col == "adj.P.Val") "FDR" else "P"

p_volc <- EnhancedVolcano(tt,
                          lab = tt$SYMBOL,
                          x = "logFC",
                          y = sig_col,
                          xlab = bquote(~Log[2] ~ "fold change"),
                          ylab = bquote(~-Log[10] ~ .(sig_label)),
                          title = cfg$title,
                          subtitle = "ACT vs Sham",
                          caption = paste("Total proteins:", nrow(tt)),
                          pCutoff = sig_cutoff,
                          FCcutoff = fc_cutoff,
                          pointSize = 1.0,
                          labSize = 3,
                          colCustom = keyvals,
                          colAlpha = 0.5,
                          legendLabSize = 7,
                          legendIconSize = 7,
                          drawConnectors = TRUE,
                          widthConnectors = 0.3,
                          max.overlaps = 15)
save_plot(p_volc, "volcano.pdf", width = 8, height = 7)
# </figure>

# <figure name="heatmap_DE">
## Heatmap of the differentially abundant proteins (batch-corrected, row-scaled)
top_ids <- c(tt_up$Protein, tt_down$Protein)
hm <- log_bc[top_ids, ]
rownames(hm) <- make.unique(tt$SYMBOL[match(top_ids, tt$Protein)])
save_base("heatmap_DE.pdf", width = 6, height = max(5, 0.15 * length(top_ids) + 2), {
  print(pheatmap(hm, scale = "row", cluster_cols = FALSE, silent = TRUE,
                 annotation_col = metadata[, c("group", "batch")],
                 annotation_colors = list(group = group_colors),
                 main = paste0("DE proteins (", sig_col, " < ", sig_cutoff, "), batch-corrected")))
})
# </figure>

# <figure name="pvalue_histogram">
## P value / FDR histogram
save_base("pvalue_histogram.pdf", width = 8, height = 4, {
  par(mfrow = c(1, 2))
  hist(tt$P.Value, breaks = 50, col = "grey", main = "P value", xlab = "P")
  hist(tt$adj.P.Val, breaks = 50, col = "grey", main = "FDR (BH)", xlab = "adj. P")
})
# </figure>

# <figure name="qq_t">
## QQ plot of moderated t
save_base("qq_t.pdf", width = 5, height = 5, {
  qqt(fit2$t[, 1], df = fit2$df.total, main = "Moderated t vs theoretical")
  abline(0, 1, col = "red")
})
# </figure>

# <figure name="md_plot">
## MD plot (logFC vs mean log intensity)
save_base("md_plot.pdf", width = 6, height = 5, {
  plotMD(fit2, coef = 1, status = tt$status[match(rownames(fit2), tt$Protein)],
         main = "ACT vs Sham", values = c("Up", "Down"), col = c("#CC0000", "#3399FF"))
  abline(h = 0, col = "grey")
})
# </figure>

# <figure name="sa_plot">
## Mean-variance trend
save_base("sa_plot.pdf", width = 6, height = 5,
          plotSA(fit2, main = "Residual SD vs mean (trend = TRUE)"))
# </figure>

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
enr <- cfg$enrichment
set.seed(cfg$seed)
gse_go <- gseGO(geneList = gene_list, ont = "ALL", keyType = "ENTREZID",
                OrgDb = org.Mm.eg.db, minGSSize = enr$min_gs_size, maxGSSize = enr$max_gs_size,
                pvalueCutoff = enr$pvalue_cutoff, pAdjustMethod = "BH", eps = 0, verbose = FALSE)
gse_go <- setReadable(gse_go, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
save(as.data.frame(gse_go), "gsea_GO")
cat("GSEA GO terms (FDR < ", enr$pvalue_cutoff, "): ", nrow(gse_go), "\n", sep = "")

if (nrow(gse_go) > 0) {
  # <figure name="gsea_GO_dotplot">
  save_plot(dotplot(gse_go, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign),
            "gsea_GO_dotplot.pdf", width = 11, height = 8)
  # </figure>
  # <figure name="gsea_GO_ridgeplot">
  save_plot(ridgeplot(gse_go, showCategory = 15) + labs(x = "signed -log10(P)"),
            "gsea_GO_ridgeplot.pdf", width = 9, height = 8)
  # </figure>
  # <figure name="gsea_GO_emapplot">
  gse_go_sim <- pairwise_termsim(gse_go)
  save_plot(emapplot(gse_go_sim, showCategory = 30), "gsea_GO_emapplot.pdf", width = 10, height = 9)
  # </figure>
  # <figure name="gsea_GO_top1">
  save_plot(gseaplot2(gse_go, geneSetID = 1, title = gse_go$Description[1], pvalue_table = TRUE),
            "gsea_GO_top1.pdf", width = 8, height = 6)
  # </figure>
}

## GSEA on KEGG (requires internet access to KEGG; can be disabled in the config)
gse_kegg <- if (!isTRUE(enr$run_kegg)) NULL else tryCatch(
  gseKEGG(geneList = gene_list, organism = "mmu", keyType = "ncbi-geneid",
          minGSSize = enr$min_gs_size, maxGSSize = enr$max_gs_size,
          pvalueCutoff = enr$pvalue_cutoff, eps = 0, verbose = FALSE),
  error = function(e) { message("KEGG skipped: ", conditionMessage(e)); NULL })
if (!is.null(gse_kegg)) {
  gse_kegg <- setReadable(gse_kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  save(as.data.frame(gse_kegg), "gsea_KEGG")
  cat("GSEA KEGG pathways (FDR < ", enr$pvalue_cutoff, "): ", nrow(gse_kegg), "\n", sep = "")
  if (nrow(gse_kegg) > 0) {
    # <figure name="gsea_KEGG_dotplot">
    save_plot(dotplot(gse_kegg, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign),
              "gsea_KEGG_dotplot.pdf", width = 11, height = 8)
    # </figure>
  }
}

## ORA of up / down proteins (same significance criterion as the DE call)
ora_sets <- list(
  up   = tt_entrez$ENTREZID[tt_entrez$status == "Up"],
  down = tt_entrez$ENTREZID[tt_entrez$status == "Down"])
ora <- compareCluster(ora_sets, fun = "enrichGO", OrgDb = org.Mm.eg.db, ont = "BP",
                      universe = universe, pvalueCutoff = enr$pvalue_cutoff, readable = TRUE)
if (!is.null(ora) && nrow(as.data.frame(ora)) > 0) {
  save(as.data.frame(ora), "ora_GO_BP_up_down")
  # <figure name="ora_GO_BP_dotplot">
  save_plot(dotplot(ora, showCategory = 10), "ora_GO_BP_dotplot.pdf", width = 9, height = 8)
  # </figure>
}

######################## Run summary ########################
# Key numbers of this run, read by the HTML report (scripts/report.py)
run_summary <- data.frame(
  key = c("proteins_input", "proteins_filtered", "missing_per_sample", "pca_var_explained",
          "sig_col", "sig_cutoff", "logfc_cutoff", "n_up", "n_down", "n_robust_noimp",
          "cor_logfc_imp_noimp", "n_gsea_go", "n_gsea_kegg", "n_ora"),
  value = c(nrow(log_int), nrow(log_filt),
            paste(names(missing_per_sample), missing_per_sample, sep = "=", collapse = "; "),
            paste0("PC", 1:4, "=", var_expl[1:4], "%", collapse = "; "),
            sig_col, sig_cutoff, fc_cutoff, nrow(tt_up), nrow(tt_down), sum(tt$robust_noimp),
            round(cor(tt$logFC, tt$noimp_logFC, use = "complete.obs"), 3),
            nrow(gse_go), if (is.null(gse_kegg)) NA else nrow(gse_kegg),
            if (is.null(ora)) 0 else nrow(as.data.frame(ora))))
save(run_summary, "run_summary")
dir.create(cfg$reports_dir, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(cfg$reports_dir, "sessionInfo.txt"))

sessionInfo()
