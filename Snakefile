# Hippocampal proteome, ACT vs Sham (McElroy et al. 2020, ijms-730741).
# Final target: reports/report.html.
#   snakemake --cores 1                 # with locally installed tools
#   snakemake --cores 1 --use-conda     # with the environments in env/
# Run from the project root (proteomics/). Parameters: config/config.yaml.

configfile: "config/config.yaml"

RESULTS = config["results_dir"]
PLOTS = config["plots_dir"]
REPORTS = config["reports_dir"]
NAME = "ACT_mcelroy_limma_ijms730741"

# Outputs always written by the analysis script. Enrichment plots, KEGG and ORA tables are
# conditional (depend on results / network access) and are therefore not declared.
TABLES = ["ACT_vs_Sham_limma_all", "ACT_vs_Sham_limma_DE", "gsea_GO", "run_summary"]
FIGURES = ["density_preprocessing", "boxplot_normalized", "pca", "sensitivity_imputation_logFC",
           "volcano", "heatmap_DE", "pvalue_histogram", "qq_t", "md_plot", "sa_plot"]


rule all:
    input:
        f"{REPORTS}/report.html",


rule prepare:
    """Extract annotation + raw iBAQ (Sham, ACT) from Table S1."""
    input:
        xlsx=config["data"]["raw_xlsx"],
        script="scripts/prepare.sh",
    output:
        config["data"]["processed"],
    log:
        "logs/prepare.log",
    conda:
        "env/prepare.yaml"
    shell:
        "bash {input.script} {input.xlsx} {output} > {log} 2>&1"


rule limma_analysis:
    """Pre-processing, limma ACT vs Sham, sensitivity without imputation, plots, GSEA/ORA."""
    input:
        csv=config["data"]["processed"],
        script=f"scripts/{NAME}.R",
        config="config/config.yaml",
    output:
        expand(f"{RESULTS}/{{t}}.csv", t=TABLES),
        f"{REPORTS}/sessionInfo.txt",
        expand(f"{PLOTS}/{{f}}.{{ext}}", f=FIGURES, ext=["pdf", "png"]),
    log:
        "logs/limma_analysis.log",
    conda:
        "env/limma.yaml"
    shell:
        "Rscript {input.script} {input.config} > {log} 2>&1"


rule report:
    """Self-contained HTML report (Python, shared html_report.py) from results/ and plots/."""
    input:
        top_table=f"{RESULTS}/ACT_vs_Sham_limma_all.csv",
        go=f"{RESULTS}/gsea_GO.csv",
        summary=f"{RESULTS}/run_summary.csv",
        session=f"{REPORTS}/sessionInfo.txt",
        plots=expand(f"{PLOTS}/{{f}}.png", f=FIGURES),
        config=["config/config.yaml"],
        envs=["env/prepare.yaml", "env/limma.yaml", "env/report.yaml"],
        scripts=["scripts/report.py", "scripts/html_report.py"],
        code=[f"scripts/{NAME}.R"],
    output:
        html=f"{REPORTS}/report.html",
    params:
        title="Hippocampal proteome: ACT vs Sham",
        subtitle=config["title"],
        design="""
Hippocampus of mice after ACT chemotherapy (doxorubicin + cyclophosphamide, then paclitaxel) or sham
treatment. Label-free LC-MS/MS quantified with MaxQuant (iBAQ). Data from the supplementary file of
McElroy et al. (2020), *Int. J. Mol. Sci.* (`ijms-730741`, Table S1).
""",
        notes="""
> **Independent re-analysis.** Only the raw iBAQ values are used. The authors' normalized values, statistics
> and protein tables are not used or compared.
>
> **Chemotherapy effect only.** The study also has MnBuOE and ACT-MnBuOE groups; they are excluded on purpose.
>
> **Replicate pair as batch.** Replicates 1–2 and 3–4 differ strongly in depth (missing values) in both groups,
> so the pair is a blocking factor in the model.
""",
        pipeline="""
1. **prepare** (`scripts/prepare.sh`): extracts annotation + raw iBAQ of the Sham and ACT samples from Table S1
   into `data/processed.csv`.
2. **limma_analysis** (`scripts/ACT_mcelroy_limma_ijms730741.R`): filtering, normalization, imputation, limma,
   sensitivity analysis without imputation, QC/DE plots, GSEA (GO, KEGG) and ORA.
3. **report** (`scripts/report.py`): this page.
""",
        input_text="""
Input is Table S1 of the supplementary file: one row per MaxQuant protein group, with raw iBAQ intensities
per sample. A value of 0 means the protein was not quantified in that sample.
""",
    log:
        "logs/report.log",
    conda:
        "env/report.yaml"
    script:
        "scripts/report.py"
