# Proteomics: chemotherapy-induced changes in the mouse hippocampal proteome

Reproducible Snakemake + R workflow for the differential abundance analysis of public label-free proteomics data
on chemotherapy-induced cognitive impairment ("chemobrain"). It is an **independent re-analysis** of the raw
MaxQuant iBAQ intensities published with McElroy et al. (2020), from the supplementary table to differential
abundance, gene set enrichment and an HTML report.

**Report:** [Hippocampal proteome: ACT vs Sham](https://arlindomatias.github.io/proteomics/reports/report.html)

The companion repository [bulk-rna-seq](https://github.com/arlindomatias/bulk-rna-seq) analyzes the
hippocampal transcriptome in the same chemotherapy model (AC-T, GSE295281) with the same conventions and report
layout.

## Dataset

- **Study:** McElroy T, Brown T, Kiffer F, Wang J, Byrum SD, Oberley-Deegan RE, Allen AR. Assessing the effects
  of redox modifier MnTnBuOE-2-PyP 5+ on cognition and hippocampal physiology following doxorubicin,
  cyclophosphamide, and paclitaxel treatment. *Int J Mol Sci* 2020;21(5):1867.
  [doi:10.3390/ijms21051867](https://doi.org/10.3390/ijms21051867), PMID 32182883.
- **Samples:** hippocampus of 4-month-old female C57BL/6 mice, about four weeks after ACT chemotherapy
  (doxorubicin + cyclophosphamide, then paclitaxel) or sham treatment. n = 4 per group.
- **Input:** Table S1 of the supplementary file (`data/ijms-730741-SI.xlsx`): MaxQuant protein groups with raw
  iBAQ intensities per sample.
- **Scope:** only **ACT vs Sham** (effect of chemotherapy). The study's MnBuOE groups are excluded on purpose, and
  the authors' normalized values and statistics are not used.

## Workflow

```
data/ijms-730741-SI.xlsx ──► prepare ──► data/processed.csv ──► limma_analysis ──► results/ + plots/ ──► report ──► reports/report.html
```

| Rule | Script | What it does |
|---|---|---|
| `prepare` | `scripts/prepare.sh` | Extracts annotation + raw iBAQ of the Sham and ACT samples from Table S1 |
| `limma_analysis` | `scripts/ACT_mcelroy_limma_ijms730741.R` | Filtering, normalization, imputation, limma, sensitivity analysis, QC/DE plots, GSEA and ORA |
| `report` | `scripts/report.py` | Self-contained HTML report (built with `scripts/html_report.py`, the same builder as bulk-rna-seq) |

### Methods

- **Missing values:** iBAQ = 0 is treated as missing. Proteins need ≥ 3 valid values in at least one group.
- **Normalization:** log2, then median centering per sample.
- **Imputation:** left-censored, MinProb-like draws from N(mean − 1.8 SD, 0.3 SD) per sample, with a fixed seed.
- **Batch:** replicates 1–2 and 3–4 differ strongly in depth in both groups. The replicate pair is a blocking
  factor in the model (`~ 0 + group + batch`); PCA and heatmap show the batch-corrected matrix
  (`limma::removeBatchEffect`, visualization only).
- **Differential abundance:** `limma` with `eBayes(trend = TRUE)`, contrast ACT − Sham (positive log2FC = higher in
  ACT). Calls use `FDR < 0.05` and `|log2FC| > 0.5`; when no protein passes FDR, the unadjusted `P < 0.01` is used
  instead and the calls are reported as exploratory.
- **Sensitivity analysis:** the same model without imputation, to flag calls that depend on imputed values.
- **Enrichment:** GSEA on GO (BP, CC, MF) and KEGG (`mmu`) with clusterProfiler, proteins ranked by
  sign(log2FC) × −log10(P); ORA of the up/down sets with the identified proteome as background.

All parameters (paths, seed, filter, imputation, cutoffs, enrichment settings, group colors) are in
[`config/config.yaml`](config/config.yaml).

### Current results

- 5,514 protein groups, 5,200 tested after filtering.
- No protein passes FDR < 0.05. With the exploratory cutoff (P < 0.01, |log2FC| > 0.5): **29 up and 9 down** in
  ACT; 20 of these 38 calls hold without imputation.
- GSEA: **104 GO terms and 12 KEGG pathways** (padj < 0.05). Enriched in ACT: mitochondrial inner membrane and
  oxidative phosphorylation, ribosome / synaptic translation, G-protein activity, and neurodegeneration pathways.
  Lower in ACT: chromatin organization and binding.

See the [report](https://arlindomatias.github.io/proteomics/reports/report.html) for figures, tables and caveats.

## Repository layout

| Folder | Contents |
|---|---|
| `config/` | `config.yaml` with all parameters |
| `data/` | Input supplementary file and the prepared table (`processed.csv`) |
| `scripts/` | All code: preparation (Bash), analysis (R), report (Python) |
| `results/` | Result tables only (`.csv`) |
| `plots/` | All figures (PDF for print, PNG embedded in the report) |
| `reports/` | `report.html` (self-contained) and `sessionInfo.txt` |
| `env/` | Conda environments (`<name>.yaml` creates the environment `<name>`) and their `<name>.linux-64.pin.txt` locks |

## Running

Requirements: Linux (x86-64), Conda and Snakemake. Each rule creates its own environment with `--use-conda`
(`prepare`, `limma`, `report`); exact versions are in `env/*.linux-64.pin.txt` (R 4.5.3, limma 3.68.5,
clusterProfiler 4.18.4). KEGG GSEA needs internet access to `rest.kegg.jp` (it can be disabled in the config).

```bash
git clone git@github.com:arlindomatias/proteomics.git
cd proteomics
snakemake -n                                             # dry run
snakemake --cores 1 --use-conda --conda-frontend conda   # full run (drop --conda-frontend if mamba is installed)
```

The R script also runs on its own: `Rscript scripts/ACT_mcelroy_limma_ijms730741.R [config/config.yaml]`.

## Data and citation

The proteomics data come from the open-access supplementary material of McElroy et al. (2020),
[doi:10.3390/ijms21051867](https://doi.org/10.3390/ijms21051867). Please cite the original authors when using
these data.

## License

[MIT](LICENSE) © 2026 Arlindo César Matias Pereira. The license covers the code and the analysis outputs of this
repository; the input data (`data/ijms-730741-SI.xlsx`) belong to the original authors and keep the license of
the source article.
