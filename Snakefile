##############################################################################
#
#   Snakemake pipeline:
#   RNA-Seq data analysis with deepTools
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 23-11-2019
#   LICENSE: GPL v3.0
#
##############################################################################

# imports
import os
import sys
import pandas as pd

# local rules
localrules: all, create_out_dir

# get path to the alignment file for a given sample
def get_input_bam(sample):
    design_table = pd.read_csv(config["design_file"], sep="\t", index_col=0)
    return design_table.at[sample,"bam"]

# get all sample names
def get_all_samples():
    design_table = pd.read_csv(config["design_file"], sep="\t", index_col=0)
    return design_table.index.values

# get the effective genome size for the genome in this analysis
def get_effective_genome_size():
    genome_sizes = \
        pd.read_csv(config["effective_genome_sizes"], sep="\t", index_col=0)
    return genome_sizes.at[config["genome_version"], "Effective size"]

##############################################################################
### Target rule with final output of the pipeline
##############################################################################

rule all:
    input:
        PNG_pca = expand(\
            os.path.join("{output_dir}", "clustering", "PCA.png"), \
                output_dir=config["output_dir"]),
        PNG_pca_transposed = expand(\
            os.path.join("{output_dir}", "clustering", "PCA_T.png"), \
                output_dir=config["output_dir"]),
        PNG_heatmap = expand(\
            os.path.join("{output_dir}", "clustering", "heatmap.png"), \
                output_dir=config["output_dir"]),
        PNG_scatterplot = expand(\
            os.path.join("{output_dir}", "clustering", "scatterplot.png"), \
                output_dir=config["output_dir"]),
        PNG_coverage_plot = expand(\
            os.path.join("{output_dir}", "clustering", "coverage_plot.png"), \
                output_dir=config["output_dir"]),
        PNG_gc_plot = expand(\
            os.path.join("{output_dir}", "GC_bias", "{sample}_gc.png"), \
                output_dir=config["output_dir"], sample=get_all_samples())

##############################################################################
### Create directories for the result
##############################################################################

rule create_out_dir:
    output:
        TMP_output = temp(os.path.join("{output_dir}", "dir_created"))
    params:
        DIR_results_dir = "{output_dir}",
        DIR_cluster_log = os.path.join("{output_dir}", "cluster_log"),
    log:
        DIR_local_log = os.path.join("{output_dir}", "local_log"),
    shell:
        """
        mkdir -p {params.DIR_results_dir}; \
        mkdir -p {params.DIR_cluster_log}; \
        mkdir -p {log.DIR_local_log}; \
        touch {output.TMP_output}
        """

##############################################################################
### Sort alignment files
##############################################################################

rule sort_bam:
    input:
        TMP_output = os.path.join("{output_dir}", "dir_created")
    output:
        BAM_sorted = \
            os.path.join("{output_dir}", "coverages", "{sample}_sorted.bam")
    params:
        BAM_path = lambda wildcards: get_input_bam(wildcards.sample),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "sort_bam_{sample}.log"),
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "sort_bam_{sample}.log")
    resources:
        threads = 4,
        mem = 20000
    benchmark:
        os.path.join("{output_dir}", "local_log", \
            "sort_bam_{sample}_benchmark.log")
    conda:
        "env_yaml/samtools.yaml"
    shell:
        """
        samtools sort -@ {resources.threads} {params.BAM_path} \
        1> {output.BAM_sorted} \
        2> {log.LOG_local_log};
        """

##############################################################################
### Index alignment file
##############################################################################

rule samtools_index:
    input:
        BAM_sorted = \
            os.path.join("{output_dir}", "coverages", "{sample}_sorted.bam")
    output:
        BAI_index = os.path.join(\
            "{output_dir}", "coverages", "{sample}_sorted.bam.bai")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "samtools_index_{sample}.log"),
        queue = "30min",
        time = "0:30:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "samtools_index_{sample}.log")
    resources:
        threads = 4,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", "local_log", \
            "samtools_index_{sample}_benchmark.log")
    conda:
        "env_yaml/samtools.yaml"
    shell:
        """
        samtools index -@ {resources.threads} {input.BAM_sorted} \
        1> {output.BAI_index} \
        2> {log.LOG_local_log};
        """

##############################################################################
### deepTools: create bigWig files from sorted and indexed bam
##############################################################################

rule bam2bigWig:
    input:
        BAM_sorted = \
            os.path.join("{output_dir}", "coverages", "{sample}_sorted.bam"),
        BAI_index = os.path.join(\
            "{output_dir}", "coverages", "{sample}_sorted.bam.bai")
    output:
        BW_sample = os.path.join(\
            "{output_dir}", "bigWig_coverages", "{sample}.bw")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "bam2bigWig_{sample}.log"),
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "bam2bigWig_{sample}.log")
    resources:
        threads = 4,
        mem = 20000
    benchmark:
        os.path.join("{output_dir}", "local_log", \
            "bam2bigWig_{sample}_benchmark.log")
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        bamCoverage \
        --bam {input.BAM_sorted} \
        --outFileName {output.BW_sample} \
        --binSize 1 \
        --outFileFormat bigwig \
        --numberOfProcessors {resources.threads} \
        2> {log.LOG_local_log};
        """

##############################################################################
### deepTools: create summary from bigWig
##############################################################################

rule bigWig_summary:
    input:
        BW_sample = expand(os.path.join(\
            "{output_dir}", "bigWig_coverages", "{sample}.bw"), \
            output_dir=config["output_dir"], sample=get_all_samples())
    output:
        NPZ_summary = os.path.join("{output_dir}", "bigWig_summary.npz")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "bigWig_summary.log"),
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "bigWig_summary.log")
    resources:
        threads = 4,
        mem = 20000
    benchmark:
        os.path.join("{output_dir}", "local_log", \
            "bigWig_summary_benchmark.log")
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        multiBigwigSummary bins \
        --bwfiles {input.BW_sample} \
        --outFileName {output.NPZ_summary} \
        --numberOfProcessors {resources.threads} \
        2> {log.LOG_local_log};
        """

##############################################################################
### deepTools: plot PCA
##############################################################################

rule plot_pca:
    input:
        NPZ_summary = os.path.join("{output_dir}", "bigWig_summary.npz")
    output:
        PNG_pca = os.path.join(\
            "{output_dir}", "clustering", "PCA.png"),
        PNG_pca_transposed = os.path.join(\
            "{output_dir}", "clustering", "PCA_T.png")
    params:
        TSV_pca_table = os.path.join(\
            "{output_dir}", "clustering", "PCA.tsv"),
        TSV_pca_transposed_table = os.path.join(\
            "{output_dir}", "clustering", "PCA_T.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "plot_pca.log"),
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "plot_pca.log")
    resources:
        threads = 1,
        mem = 20000
    benchmark:
        os.path.join("{output_dir}", "local_log", \
            "plot_pca_benchmark.log")
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        plotPCA \
        --corData {input.NPZ_summary} \
        --plotFile {output.PNG_pca} \
        --outFileNameData {params.TSV_pca_table} \
        --ntop 1000 \
        2> {log.LOG_local_log};
        plotPCA --transpose \
        --corData {input.NPZ_summary} \
        --plotFile {output.PNG_pca_transposed} \
        --outFileNameData {params.TSV_pca_transposed_table} \
        --ntop 1000 \
        2> {log.LOG_local_log};
        """

##############################################################################
### deepTools: plot correlation
##############################################################################

rule correlation_plot:
    input:
        NPZ_summary = "{output_dir}/bigWig_summary.npz"
    output:
        PNG_heatmap = os.path.join(\
            "{output_dir}", "clustering", "heatmap.png"),
        PNG_scatterplot = os.path.join(\
            "{output_dir}", "clustering", "scatterplot.png")
    params:
        TSV_heatmap_table = os.path.join(\
            "{output_dir}", "clustering", "heatmap.tsv"),
        TSV_scatterplot_table = os.path.join(\
            "{output_dir}", "clustering", "scatterplot.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "correlation_plot.log"),
        queue = "30min",
        time = "0:30:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "correlation_plot.log")
    resources:
        threads = 1,
        mem = 20000
    benchmark:
        os.path.join("{output_dir}", "local_log", \
            "correlation_plot_benchmark.log")
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        plotCorrelation \
        --corData {input.NPZ_summary} \
        --corMethod pearson \
        --whatToPlot heatmap \
        --plotFile {output.PNG_heatmap} \
        --outFileCorMatrix {params.TSV_heatmap_table} \
        2> {log.LOG_local_log};
        plotCorrelation \
        --corData {input.NPZ_summary} \
        --corMethod pearson \
        --whatToPlot scatterplot \
        --plotFile {output.PNG_scatterplot} \
        --outFileCorMatrix {params.TSV_scatterplot_table} \
        2> {log.LOG_local_log};
        """

##############################################################################
### deepTools: plot coverage information
##############################################################################

rule plot_coverage:
    input:
        BAM_sorted = expand(os.path.join(\
            "{output_dir}", "coverages", "{sample}_sorted.bam"), \
            output_dir=config["output_dir"], sample=get_all_samples()),
        BAI_index = expand(os.path.join(\
            "{output_dir}", "coverages", "{sample}_sorted.bam.bai"), \
            output_dir=config["output_dir"], sample=get_all_samples())
    output:
        PNG_coverage_plot = \
            os.path.join("{output_dir}", "clustering", "coverage_plot.png")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "plot_coverage.log"),
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "plot_coverage.log")
    resources:
        threads = 4,
        mem = 50000
    benchmark:
        os.path.join("{output_dir}", "local_log", \
            "plot_coverage_benchmark.log")
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        plotCoverage \
        --bamfiles {input.BAM_sorted} \
        --plotFile {output.PNG_coverage_plot} \
        --numberOfProcessors {resources.threads} \
        2> {log.LOG_local_log};
        """

##############################################################################
### convert fasta genome to .2bit format
##############################################################################

rule fa_to_2bit:
    output:
        TWOBIT_genome_2bit = os.path.join("{output_dir}", "genome.2bit")
    params:
        FASTA_genome = config["genome"],
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "fa_to_2bit.log"),
        queue = "30min",
        time = "0:30:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "fa_to_2bit.log")
    resources:
        threads = 1,
        mem = 20000
    benchmark:
        os.path.join("{output_dir}", "local_log", \
            "fa_to_2bit_benchmark.log")
    conda:
        "env_yaml/ucsc.yaml"
    shell:
        """
        faToTwoBit {params.FASTA_genome} {output.TWOBIT_genome_2bit} \
        2> {log.LOG_local_log};
        """

##############################################################################
### deepTools: GC bias plots
##############################################################################

rule GC_plot:
    input:
        BAM_sorted = \
            os.path.join("{output_dir}", "coverages", "{sample}_sorted.bam"),
        BAI_index = os.path.join(\
            "{output_dir}", "coverages", "{sample}_sorted.bam.bai"),
        TWOBIT_genome_2bit = os.path.join("{output_dir}", "genome.2bit")
    output:
        PNG_gc_plot = \
            os.path.join("{output_dir}", "GC_bias", "{sample}_gc.png")
    params:
        TSV_gc_tsv = \
            os.path.join("{output_dir}", "GC_bias", "{sample}_gc.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "GC_plot_{sample}.log"),
        effective_genome_size = get_effective_genome_size(),
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "GC_plot_{sample}.log")
    resources:
        threads = 4,
        mem = 20000
    benchmark:
        os.path.join("{output_dir}", "local_log", \
            "GC_plot_{sample}_benchmark.log")
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        computeGCBias \
        --bamfile {input.BAM_sorted} \
        --effectiveGenomeSize {params.effective_genome_size} \
        --genome {input.TWOBIT_genome_2bit} \
        --numberOfProcessors {resources.threads} \
        --GCbiasFrequenciesFile {params.TSV_gc_tsv} \
        --biasPlot {output.PNG_gc_plot} \
        --plotFileFormat png \
        2> {log.LOG_local_log};
        """
