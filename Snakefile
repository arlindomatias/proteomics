# Hippocampal proteome, ACT vs Sham (McElroy et al. 2020, ijms-730741).
#   snakemake --cores 1                 # with locally installed tools (final target: HTML report)
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
        f"{REPORTS}/{NAME}.html",


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
        f"{RESULTS}/session_info.txt",
        expand(f"{PLOTS}/{{f}}.{{ext}}", f=FIGURES, ext=["pdf", "png"]),
    log:
        "logs/limma_analysis.log",
    conda:
        "env/r_analysis.yaml"
    shell:
        "Rscript {input.script} {input.config} > {log} 2>&1"


rule report:
    """Self-contained HTML report built from the result tables and PNG figures."""
    input:
        rmd=f"scripts/report_{NAME}.Rmd",
        config="config/config.yaml",
        tables=expand(f"{RESULTS}/{{t}}.csv", t=TABLES),
        session=f"{RESULTS}/session_info.txt",
        figures=expand(f"{PLOTS}/{{f}}.png", f=FIGURES),
    output:
        f"{REPORTS}/{NAME}.html",
    log:
        "logs/report.log",
    conda:
        "env/report.yaml"
    shell:
        """
        Rscript -e 'rmarkdown::render("{input.rmd}", output_file = basename("{output}"),
            output_dir = dirname("{output}"), knit_root_dir = getwd(),
            intermediates_dir = tempdir(), params = list(config = "{input.config}"))' > {log} 2>&1
        """
