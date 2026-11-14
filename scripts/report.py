"""HTML report for the limma proteomics project: input → QC → DE → sensitivity → GSEA/ORA.

Run through Snakemake (`script:` directive, rule `report`). Project-specific text comes from
`snakemake.params` and cutoffs from the config. It only reads results/, plots/ and reports/sessionInfo.txt;
every number in the text is computed from those tables.
"""
import os
import sys

import pandas as pd

# Shared report module: scripts/html_report.py (identical copy of bulk-rna-seq/scripts/html_report.py).
sys.path.insert(0, snakemake.scriptdir)  # noqa: F821
from html_report import Report  # noqa: E402

inp, out, prm, cfg = snakemake.input, snakemake.output, snakemake.params, snakemake.config  # noqa: F821
P = cfg["plots_dir"]
ref, trt = cfg["groups"]["reference"], cfg["groups"]["treatment"]
de_cfg, pre, enr = cfg["de"], cfg["preprocessing"], cfg["enrichment"]


def read_csv(path):
    try:
        return pd.read_csv(path)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


tt = read_csv(inp.top_table)
go = read_csv(inp.go)
kegg = read_csv(os.path.join(cfg["results_dir"], "gsea_KEGG.csv"))      # optional (needs KEGG access)
ora = read_csv(os.path.join(cfg["results_dir"], "ora_GO_BP_up_down.csv"))  # optional (may be empty)
summ = read_csv(inp.summary).set_index("key")["value"]

sig_col = summ["sig_col"]
sig_cut = float(summ["sig_cutoff"])
lfc_cut = float(summ["logfc_cutoff"])
sig_name = "FDR" if sig_col == "adj.P.Val" else "unadjusted P"
used_fallback = sig_col == "P.Value"

up = tt[tt["status"] == "Up"]
down = tt[tt["status"] == "Down"]
deg = pd.concat([up, down])
robust = deg[deg["robust_noimp"]]
fdr_only = int((tt["adj.P.Val"] < de_cfg["fdr_cutoff"]).sum())

missing = dict(kv.split("=") for kv in summ["missing_per_sample"].split("; "))
samples = pd.DataFrame({"sample": list(missing), "missing values": [int(v) for v in missing.values()]})
samples["group"] = samples["sample"].str.split("_").str[0]
samples["replicate"] = samples["sample"].str.split("_").str[1].astype(int)
samples["batch (replicate pair)"] = ["B1" if r <= 2 else "B2" for r in samples["replicate"]]
samples["quantified"] = int(summ["proteins_input"]) - samples["missing values"]
n_ref, n_trt = (samples["group"] == ref).sum(), (samples["group"] == trt).sum()


def sig(df, col="p.adjust"):
    return df[df[col] < enr["pvalue_cutoff"]] if len(df) and col in df else df.iloc[0:0]


go_sig, kegg_sig = sig(go), sig(kegg)

r = Report(prm.title, subtitle=prm.subtitle, code_files=list(inp.code))
r.kpis([
    (len(samples), "samples"),
    (f"{int(summ['proteins_filtered']):,}", "proteins tested"),
    (f"{len(up):,}", f"up in {trt}"),
    (f"{len(down):,}", f"down in {trt}"),
    (f"{len(robust):,} / {len(deg):,}", "DE calls robust without imputation"),
    (f"{len(go_sig):,} / {len(kegg_sig):,}", f"GSEA terms GO / KEGG (padj < {enr['pvalue_cutoff']})"),
])

# ---------------------------------------------------------------- design
r.h2("Study design")
r.md(prm.design)
r.md(f"* **{ref}**: {n_ref} samples.\n* **{trt}**: {n_trt} samples.")
if prm.notes:
    r.md(prm.notes)

r.h2("Pipeline")
r.md(prm.pipeline)
r.md(f"""
Shared settings:

* **Filter:** proteins with ≥ {pre['min_valid_per_group']} valid (non-zero) iBAQ values in at least one group.
* **Normalization:** log2, then median centering per sample.
* **Imputation:** left-censored (MinProb-like) draws from N(mean − {pre['impute_shift']} SD, {pre['impute_width']} SD)
  per sample, seed {cfg['seed']}.
* **limma:** `~ 0 + group + batch` (batch = replicate pair 1–2 vs 3–4), `eBayes(trend = TRUE)`, contrast
  {trt} − {ref}. Positive log2FC = higher in **{trt}**.
* **DE call:** `FDR < {de_cfg['fdr_cutoff']}` and `|log2FC| > {lfc_cut}`; if no protein passes FDR, the
  unadjusted `P < {de_cfg['nominal_cutoff']}` is used instead (same fold-change cutoff).
* **GSEA:** proteins ranked by sign(log2FC) × −log10(P); GO (BP, CC, MF) and KEGG (`mmu`) with clusterProfiler,
  gene sets {enr['min_gs_size']}–{enr['max_gs_size']}, Benjamini-Hochberg, fixed seed. Positive NES = enriched in {trt}.
* **ORA:** GO BP of the up and down sets against the identified proteome (not the whole genome).
""")

# ---------------------------------------------------------------- input
r.h2("Input data")
r.md(prm.input_text)
r.md(f"""
**{int(summ['proteins_input']):,}** protein groups are in the table; **{int(summ['proteins_filtered']):,}** pass the
valid-value filter. Replicates 1–2 have many more missing values than replicates 3–4 in both groups.
""")
r.table(samples[["sample", "group", "batch (replicate pair)", "quantified", "missing values"]],
        "Proteins quantified per sample (iBAQ > 0).")

# ---------------------------------------------------------------- QC
r.h2("Quality control")
r.md("""
These plots check the data before testing. After median centering the samples should have similar
distributions, and replicates of the same group should be close in the PCA.
""")
r.h3("Intensity distributions")
r.figure(f"{P}/density_preprocessing.png",
         "Density of log2 iBAQ per sample after filtering, after normalization and after imputation. "
         "The imputed values form the low-intensity shoulder.", optional=True)
r.figure(f"{P}/boxplot_normalized.png", "Normalized log2 iBAQ per sample (before imputation), colored by group.",
         optional=True)
r.h3("Sample similarity")
pcs = summ["pca_var_explained"].replace("; ", ", ")
r.md(f"""
Replicates 1–2 and 3–4 form two batches (replicate pairs) with very different depth, each with 2 {ref} and
2 {trt} samples. Without correction the batch is the largest source of variation in the PCA. The limma model
estimates it as a term (`~ 0 + group + batch`); for the plots below the batch is removed with
`limma::removeBatchEffect` (group effect protected), so they show the variation left after the batch.

Variance explained after batch correction: {pcs}.
""")
r.figure(f"{P}/pca.png", "PCA of the imputed log2 intensities after removing the batch. Color = group, "
         "shape = batch (replicate pair).", width="80%", optional=True)

# ---------------------------------------------------------------- DE
r.h2("Differential abundance")
crit = f"`{sig_name} < {sig_cut}` and `|log2FC| > {lfc_cut}`"
fallback_txt = (f"No protein passes `FDR < {de_cfg['fdr_cutoff']}` ({fdr_only} without the fold-change cutoff), "
                f"so the unadjusted `P < {sig_cut}` is used. These calls are **exploratory**: with ~{len(tt):,} "
                f"tests, about {len(tt) * sig_cut:.0f} proteins are expected below `P < {sig_cut}` by chance alone."
                if used_fallback else
                f"{fdr_only:,} proteins have `FDR < {de_cfg['fdr_cutoff']}` without the fold-change cutoff.")
r.md(f"""
{len(tt):,} proteins were tested. **{len(deg):,} proteins** pass {crit}: **{len(up):,} up** and
**{len(down):,} down** in {trt}.

{fallback_txt}
""")
r.figure(f"{P}/volcano.png", f"Volcano plot: log2 fold change vs {sig_name}. Red = up in {trt}, blue = down, "
         "gray = not significant.", optional=True)
r.figure(f"{P}/md_plot.png", "MD plot: mean log2 intensity vs log2 fold change (red = up, blue = down).",
         width="75%", optional=True)
r.figure(f"{P}/heatmap_DE.png", "Differentially abundant proteins: imputed log2 intensities after removing the "
         "batch, scaled per row.",
         width="70%", optional=True)

r.h3("P-value diagnostics")
r.md("A well-behaved test has a flat p-value histogram with a peak near zero, and a QQ plot that follows the "
     "diagonal for most proteins. `plotSA` shows the intensity-dependent variance trend used by `eBayes(trend = TRUE)`.")
r.figure(f"{P}/pvalue_histogram.png", "Histograms of p-values and BH-adjusted p-values.", width="80%", optional=True)
r.figure(f"{P}/qq_t.png", "QQ plot of the moderated t statistics.", width="55%", optional=True)
r.figure(f"{P}/sa_plot.png", "Residual standard deviation vs mean intensity, with the fitted trend.",
         width="65%", optional=True)

r.h3("Differentially abundant proteins")
cols = ["SYMBOL", "protein_name", "logFC", "AveExpr", "P.Value", "adj.P.Val",
        "n_missing_Sham", "n_missing_ACT", "robust_noimp"]
fmt = lambda df: df.sort_values("P.Value")[cols].assign(  # noqa: E731
    protein_name=lambda d: d["protein_name"].astype(str).str.slice(0, 50),
    robust_noimp=lambda d: d["robust_noimp"].map({True: "yes", False: "no"}))
r.table(fmt(up), f"Proteins up in {trt} (by P). `n_missing_*` = replicates without a value (imputed); "
        "`robust_noimp` = the call holds without imputation.")
r.table(fmt(down), f"Proteins down in {trt} (by P).")
r.md("The full table is `results/ACT_vs_Sham_limma_all.csv`.")

# ---------------------------------------------------------------- sensitivity
r.h2("Sensitivity to imputation")
dep = deg[~deg["robust_noimp"]]
absent = dep[dep["noimp_logFC"].isna()]
r.md(f"""
The same limma model was refit on the normalized matrix **without imputation** (limma skips missing values;
proteins absent from a whole group get no estimate). The log2FC with and without imputation correlate with
r = {summ['cor_logfc_imp_noimp']}.

**{len(robust):,} of {len(deg):,}** DE calls hold without imputation (same direction, same cutoffs). Of the
{len(dep):,} that do not, {len(absent)} have no estimate without imputation (missing in a whole group): they are
presence/absence patterns rather than quantitative changes, and their large fold changes come mostly from the
imputed values.
""")
r.figure(f"{P}/sensitivity_imputation_logFC.png",
         "log2FC with vs without imputation. Red = DE call robust without imputation.", width="80%", optional=True)
if len(dep):
    r.table(dep.sort_values("P.Value")[["SYMBOL", "status", "logFC", "P.Value", "noimp_logFC", "noimp_P.Value",
                                        "n_missing_Sham", "n_missing_ACT"]],
            "DE calls that do not hold without imputation.")

# ---------------------------------------------------------------- GSEA
r.h2("Gene set enrichment (GSEA)")
r.md("""
GSEA uses all tested proteins, so it also catches pathways where many proteins change a little. This is
where the clearest signal of this dataset is, since few single proteins pass the cutoffs.
""")


def top_by_direction(df, cols, n=10):
    if not len(df):
        return df
    pos = df[df["NES"] > 0].sort_values("p.adjust").head(n)
    neg = df[df["NES"] < 0].sort_values("p.adjust").head(n)
    return pd.concat([pos, neg])[cols]


r.h3("GO")
if len(go):
    n_up = int((go_sig["NES"] > 0).sum())
    r.md(f"**{len(go_sig):,} GO terms** have `padj < {enr['pvalue_cutoff']}` ({n_up} enriched in {trt}, "
         f"{len(go_sig) - n_up} in {ref}).")
r.figure(f"{P}/gsea_GO_dotplot.png", "Top GO terms, split by direction (activated: NES > 0, higher in ACT; "
         "suppressed: NES < 0).", optional=True)
r.figure(f"{P}/gsea_GO_ridgeplot.png", "Distribution of the ranking statistic of the core proteins of each term.",
         optional=True)
r.figure(f"{P}/gsea_GO_emapplot.png", "Enrichment map: terms that share proteins are linked.", optional=True)
r.figure(f"{P}/gsea_GO_top1.png", "Running enrichment score of the most significant GO term.", optional=True)
if len(go_sig):
    r.table(top_by_direction(go_sig, [c for c in ["ONTOLOGY", "ID", "Description", "setSize", "NES", "p.adjust"]
                                      if c in go_sig]),
            f"Top 10 GO terms per direction (NES > 0: enriched in {trt}).")

r.h3("KEGG")
if len(kegg):
    n_up = int((kegg_sig["NES"] > 0).sum())
    r.md(f"**{len(kegg_sig):,} KEGG pathways** have `padj < {enr['pvalue_cutoff']}` ({n_up} enriched in {trt}, "
         f"{len(kegg_sig) - n_up} in {ref}).")
else:
    r.md("KEGG GSEA was not run (disabled in the config or KEGG not reachable).")
r.figure(f"{P}/gsea_KEGG_dotplot.png", "Top KEGG pathways, split by direction.", optional=True)
if len(kegg_sig):
    r.table(top_by_direction(kegg_sig, [c for c in ["ID", "Description", "setSize", "NES", "p.adjust"] if c in kegg_sig]),
            f"Top 10 KEGG pathways per direction (NES > 0: enriched in {trt}).")

r.h2("Over-representation (ORA)")
if len(ora):
    r.md(f"GO BP terms over-represented in the up and down sets (padj < {enr['pvalue_cutoff']}), with the "
         "identified proteome as background.")
    r.figure(f"{P}/ora_GO_BP_dotplot.png", "ORA of the up and down sets.", optional=True)
    r.table(ora.sort_values("p.adjust")[["Cluster", "Description", "GeneRatio", "p.adjust", "geneID"]],
            "Significant GO BP terms.", max_rows=30)
else:
    r.md(f"No GO BP term is over-represented in the up ({len(up)}) or down ({len(down)}) sets at "
         f"`padj < {enr['pvalue_cutoff']}`: the sets are small compared with the ~{len(tt):,}-protein background.")

# ---------------------------------------------------------------- appendix
r.appendix(
    outputs=[
        ("data/processed.csv", "Annotation + raw iBAQ of the Sham and ACT samples (from Table S1)"),
        ("results/ACT_vs_Sham_limma_all.csv", "limma results for every tested protein (with and without imputation)"),
        ("results/ACT_vs_Sham_limma_DE.csv", "Differentially abundant proteins"),
        ("results/gsea_GO.csv", "GSEA, GO (BP, CC, MF)"),
        ("results/gsea_KEGG.csv", "GSEA, KEGG"),
        ("results/run_summary.csv", "Key numbers of the run (read by this report)"),
    ],
    env_files=sorted(inp.envs),
    config_files=list(inp.config),
    session_file=inp.session,
)

r.write(out.html)
