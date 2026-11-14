# proteomics

Differential protein abundance analysis of public label-free proteomics datasets (mouse). The folder follows
the same organization as `../rna_seq`: **one self-contained R script per dataset**, run interactively section by
section, with no shared package or pipeline framework.

## Layout

```
proteomics/
├── CLAUDE.md
├── ijms-730741-SI.xlsx        # raw supplementary file (input, never edit)
├── prepare.sh                 # extracts the analysis matrix from the xlsx -> processed.csv
├── processed.csv              # tidy input for R (generated; rerun prepare.sh to rebuild)
└── <model>_<author>_<method>_<accession>.R   # one analysis script per dataset
```

- Script naming mirrors `rna_seq` (e.g. `ACT_cronin_DESeq2_GSE295281.R`):
  `<model>_<firstauthor>_<method>_<accession>.R`; for the current dataset: `ACT_mcelroy_limma_ijms730741.R`.
  Use the PRIDE accession (`PXD...`) when one exists, otherwise the article id (`ijms-730741`).
- Data preparation (xlsx/raw -> CSV) is done in Bash (`prepare.sh`), with inline Python only when Bash is not
  enough. R starts from the prepared CSV. Language priority: Bash > Python > R.
- Each new dataset gets its own `prepare_<accession>.sh` if more than one dataset lives here.
- Supplementary/raw files are inputs: never modify them in place.

## Current dataset: `ijms-730741-SI.xlsx`

- Source: McElroy et al. (2020), "Assessing the Effects of Redox Modifier MnTnBuOE-2-PyP 5+ on Cognition and
  Hippocampal Physiology Following Doxorubicin, Cyclophosphamide, and Paclitaxel Treatment", *Int. J. Mol. Sci.*
  (published 9 March 2020). Supplementary file id: `ijms-730741`.
- Mouse **hippocampus**; ACT = doxorubicin (A) + cyclophosphamide (C) + paclitaxel (T) chemotherapy, the same
  model as `../rna_seq/ACT_cronin_DESeq2_GSE295281.R`.
- MaxQuant output, **iBAQ** intensities. Groups in the study: Sham (sample ids 100-103),
  ACT (200-203), MnBuOE (300-303), ACT-MnBuOE (400-403); n = 4 per group.
- Table S1: annotation + raw iBAQ + normalized values, followed by two limma result blocks (not row-aligned).
  Only the raw iBAQ is used. Tables S2/S3 (proteins from the article's Tables 1/2) are not used.
- **Independent analysis.** Do not use, extract or compare against the authors' normalized values, limma
  results or protein tables; everything is derived from the raw iBAQ.
- `prepare.sh` keeps only annotation + raw iBAQ for Sham and ACT and renames samples to `Sham_1..4`, `ACT_1..4`.
  `processed.csv` has 5514 proteins; columns: `Protein, Protein names, Gene names, Fasta headers,
  Mol. weight [kDa], Sequence length, Sham_1..Sham_4, ACT_1..ACT_4`.
- Missing values are encoded as `0` in iBAQ: convert to `NA` before log transformation.
- `Gene names` may hold several symbols separated by `;` (protein groups): keep the first for annotation
  and ID mapping, or collapse duplicates explicitly.
- Main contrast: **ACT vs Sham** (positive logFC = higher in ACT).
- **Scope: chemotherapy effect only.** MnBuOE and ACT-MnBuOE are excluded on purpose; do not add them back
  to `processed.csv` or to the design/contrasts.

## Analysis script structure

Same section order as the `rna_seq` limma scripts, with `######################## Section ########################`
headers and `## Subsection` for plots:

1. **Packages** — `tidyverse`, `limma`, `clusterProfiler`, `org.Mm.eg.db`, `enrichplot`, `EnhancedVolcano`,
   `pheatmap`, `DOSE`, `pathview` (+ `DEP`/`imputeLCMD`/`vsn` if needed for proteomics steps).
2. **Adjusts** — plotting parameters, palette (`#3399FF`, `#CC0000`, `#66A61E`), helper `save()` for CSV.
3. **Import Data** — `read.csv("processed.csv", check.names = FALSE)`; build sample metadata (`group` factor,
   levels `Sham`, `ACT`) from column names.
4. **Pre-processing** (proteomics-specific)
   - remove contaminants / reverse hits if flagged; `0 -> NA`; `log2` transform;
   - filter by valid values (e.g. >= 3 of 4 in at least one group);
   - normalization (median centering or VSN) — check density/boxplots before and after;
   - imputation of the remaining missing values (MinProb/QRILC for MNAR, kNN for MAR), stated explicitly.
5. **Differential expression** — `limma`: `model.matrix(~ 0 + group)`, `makeContrasts(ACT - Sham)`,
   `lmFit` -> `contrasts.fit` -> `eBayes` (consider `trend = TRUE`), `topTable(number = Inf)`.
   Significance rule: BH FDR < 0.05 and |logFC| > 0.5; **if nothing passes FDR, use unadjusted P < 0.01**
   (same |logFC|). The script picks the criterion automatically and prints which one was used.
6. **Visualization** — volcano, heatmap of DE proteins (row-scaled), boxplot, p-value/FDR histogram,
   density, QQ plot of t, MD plot, PCA of samples.
7. **Enrichment** — map `SYMBOL -> ENTREZID` with `bitr(OrgDb = org.Mm.eg.db)`; GSEA on GO (ALL/BP/MF/CC)
   and KEGG (`organism = "mmu"`) with the ranked logFC list; ORA of up/down sets; dotplot, ridgeplot,
   cnetplot, emapplot, gseaplot2. Use the identified proteome as background (universe) for ORA, not the
   whole genome.

## Conventions

- All code, comments and section headers in **English** (translate Portuguese when copying from `rna_seq`).
- Relative paths only (no `C:/Users/...` or absolute home paths); outputs go next to the script or in
  `results/` (CSV tables) and `figures/` (plots), both git-ignored if large.
- Run R scripts with `Rscript <script>.R` for a full run; interactive use in RStudio is also expected.
