##############################################################################
#
#   Snakemake pipeline:
#   deepTools
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 23-11-2019
#   LICENSE: GPL v3.0
#
##############################################################################


'''
# imports
import sys
import os

# local rules
localrules: create_output_dir, all

##############################################################################
### Target rule with final output of the pipeline
##############################################################################

rule all:
    input:
        TXT_final_results = \
            expand(os.path.join("{output_dir}", "results.txt"),
                output_dir=config["output_dir"])

##############################################################################
### Create directories for the result
##############################################################################

rule create_output_dir:
    output:
        TMP_output = temp(os.path.join("{output_dir}", "dir_created"))
    params:
        DIR_random_samples = os.path.join("{output_dir}", "random_samples"),
        DIR_results_dir = "{output_dir}",
        DIR_cluster_log = os.path.join("{output_dir}", "cluster_log"),
    log:
        DIR_local_log = os.path.join("{output_dir}", "local_log"),
    shell:
        """
        mkdir -p {params.DIR_results_dir}; \
        mkdir -p {params.DIR_random_samples}; \
        mkdir -p {params.DIR_cluster_log}; \
        mkdir -p {log.DIR_local_log}; \
        touch {output.TMP_output}
        """

##############################################################################
### Sample some random data
##############################################################################

rule generate_files:
    input:
        TMP_output = os.path.join("{output_dir}", "dir_created"),
        SCRIPT = \
            os.path.join(config["src_dir"], "mb_random_sample.py")
    output:
        TXT_random_sample = \
            os.path.join("{output_dir}", "random_samples", "{file}")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "generate_files_{file}.log"),
        queue = "30min",
        time = "0:05:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "generate_files_{file}.log"),
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}",
            "cluster_log", "generate_files_{file}_benchmark.log")
    conda:
        "packages.yaml"
    singularity:
        ""
    shell:
        """
        python {input.SCRIPT} \
        --outfile {output.TXT_random_sample} \
        &> {log.LOG_local_log};
        """
'''







'''
Pipeline to plot information about genomic coverages.

Maciek Bak
'''

import os
import sys
import pandas as pd

localrules: final, create_out_dir

def get_input_bam(sample):
    design_table = pd.read_csv(config["design_file"],sep="\t", index_col=0)
    return design_table.at[sample,"bam"]

def get_all_samples():
    design_table = pd.read_csv(config["design_file"],sep="\t", index_col=0)
    return design_table.index.values

#################################################################################
### Target rule with final outfiles
#################################################################################

rule final:
    input:
        pca = expand("{output_dir}/clustering_info/PCA.png", output_dir=config["output_dir"]),
        pca_transposed = expand("{output_dir}/clustering_info/PCA_T.png", output_dir=config["output_dir"]),
        heatmap = expand("{output_dir}/clustering_info/heatmap.png", output_dir=config["output_dir"]),
        scatterplot = expand("{output_dir}/clustering_info/scatterplot.png", output_dir=config["output_dir"]),
        coverage_plot = expand("{output_dir}/coverage_plot.png", output_dir=config["output_dir"]),
        gc_plot = expand("{output_dir}/GC_bias/{sample}_gc.png", output_dir=config["output_dir"], sample=get_all_samples()),

#################################################################################
### Create directories for the result and copy config immediately
#################################################################################

rule create_out_dir:
    output:
        TMP_output = temp(os.path.join("{output_dir}", "dir_created"))
    params:
        main_dir = "{output_dir}",
        LOG_cluster_log = "{output_dir}/cluster_log",
    log:
        LOG_local_log = "{output_dir}/local_log"
    shell:
        """
        mkdir -p {params.main_dir}; \
        mkdir -p {params.cluster_log}; \
        mkdir -p {log.LOG_local_log}; \
        touch {output.TMP_output}
        """

#################################################################################
### Sort alignment files
#################################################################################

rule sort_bam:
    input:
        TMP_output = os.path.join("{output_dir}", "dir_created"),
        bam = lambda wildcards: get_input_bam(wildcards.sample)
    output:
        bam = "{output_dir}/bam_coverages/{sample}_sorted.bam"
    params:
        LOG_cluster_log = "{output_dir}/cluster_log/sort_bam_{sample}.log",
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = "{output_dir}/local_log/sort_bam_{sample}.log",
    resources:
        threads = 8,
        mem = 20000
    benchmark:
        "{output_dir}/cluster_log/sort_bam_{sample}_benchmark.log"
    conda:
        "env_yaml/samtools.yaml"
    shell:
        """
        samtools sort -@ {resources.threads} {input.bam} 1> {output.bam} \
        2> {log.LOG_local_log};
        """

#################################################################################
### Index alignment file
#################################################################################

rule samtools_index:
    input:
        bam = "{output_dir}/bam_coverages/{sample}_sorted.bam"
    output:
        bai = "{output_dir}/bam_coverages/{sample}_sorted.bam.bai"
    params:
        LOG_cluster_log = "{output_dir}/cluster_log/samtools_index_{sample}.log",
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = "{output_dir}/local_log/samtools_index_{sample}.log",
    resources:
        threads = 8,
        mem = 5000
    benchmark:
        "{output_dir}/cluster_log/samtools_index_{sample}_benchmark.log"
    conda:
        "env_yaml/samtools.yaml"
    shell:
        """
        samtools index -@ {resources.threads} {input.bam} 1> {output.bai} \
        2> {log.LOG_local_log};
        """

#################################################################################
### deepTools: create bigWig files from sorted bam
#################################################################################

rule bam2bigWig:
    input:
        bam = "{output_dir}/bam_coverages/{sample}_sorted.bam",
        bai = "{output_dir}/bam_coverages/{sample}_sorted.bam.bai"
    output:
        bigWig = "{output_dir}/bigWig_coverages/{sample}.bw"
    params:
        LOG_cluster_log = "{output_dir}/cluster_log/bam2bigWig_{sample}.log",
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = "{output_dir}/local_log/bam2bigWig_{sample}.log",
    resources:
        threads = 8,
        mem = 20000
    benchmark:
        "{output_dir}/cluster_log/bam2bigWig_{sample}_benchmark.log"
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        bamCoverage \
        --bam {input.bam} \
        --outFileName {output.bigWig} \
        --binSize 1 \
        --outFileFormat bigwig \
        --numberOfProcessors {resources.threads} \
        2> {log.LOG_local_log};
        """

#################################################################################
### deepTools: create summary from bigWig
#################################################################################

rule bigWig_summary:
    input:
        bigWig = expand("{output_dir}/bigWig_coverages/{sample}.bw", output_dir=config["output_dir"], sample=get_all_samples())
    output:
        summary = "{output_dir}/bigWig_summary.npz"
    params:
        LOG_cluster_log = "{output_dir}/cluster_log/bigWig_summary.log",
        queue = "6hours",
        time = "2:00:00"
    log:
        LOG_local_log = "{output_dir}/local_log/bigWig_summary.log",
    resources:
        threads = 8,
        mem = 20000
    benchmark:
        "{output_dir}/cluster_log/bigWig_summary_benchmark.log"
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        multiBigwigSummary bins \
        --bwfiles {input.bigWig} \
        --outFileName {output.summary} \
        --numberOfProcessors {resources.threads} \
        2> {log.LOG_local_log};
        """

#################################################################################
### deepTools: plot PCA
#################################################################################

rule plot_pca:
    input:
        summary = "{output_dir}/bigWig_summary.npz"
    output:
        pca = "{output_dir}/clustering_info/PCA.png",
        pca_transposed = "{output_dir}/clustering_info/PCA_T.png"
    params:
        pca_table = "{output_dir}/clustering_info/PCA.tsv",
        pca_transposed_table = "{output_dir}/clustering_info/PCA_T.tsv",
        LOG_cluster_log = "{output_dir}/cluster_log/plot_pca.log",
        queue = "6hours",
        time = "1:00:00"
    log:
        LOG_local_log = "{output_dir}/local_log/plot_pca.log",
    resources:
        threads = 1,
        mem = 20000
    benchmark:
        "{output_dir}/cluster_log/plot_pca_benchmark.log"
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        plotPCA \
        --corData {input.summary} \
        --plotFile {output.pca} \
        --outFileNameData {params.pca_table} \
        --ntop 1000 \
        2> {log.LOG_local_log};
        plotPCA --transpose \
        --corData {input.summary} \
        --plotFile {output.pca_transposed} \
        --outFileNameData {params.pca_transposed_table} \
        --ntop 1000 \
        2> {log.LOG_local_log};
        """

#################################################################################
### deepTools: plot correlation
#################################################################################

rule correlation_plot:
    input:
        summary = "{output_dir}/bigWig_summary.npz"
    output:
        heatmap = "{output_dir}/clustering_info/heatmap.png",
        scatterplot = "{output_dir}/clustering_info/scatterplot.png"
    params:
        heatmap_table = "{output_dir}/clustering_info/heatmap.tsv",
        scatterplot_table = "{output_dir}/clustering_info/scatterplot.tsv",
        LOG_cluster_log = "{output_dir}/cluster_log/correlation_plot.log",
        queue = "6hours",
        time = "1:00:00"
    log:
        LOG_local_log = "{output_dir}/local_log/correlation_plot.log",
    resources:
        threads = 1,
        mem = 20000
    benchmark:
        "{output_dir}/cluster_log/correlation_plot_benchmark.log"
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        plotCorrelation \
        --corData {input.summary} \
        --corMethod pearson \
        --whatToPlot heatmap \
        --plotFile {output.heatmap} \
        --outFileCorMatrix {params.heatmap_table} \
        2> {log.LOG_local_log};
        plotCorrelation \
        --corData {input.summary} \
        --corMethod pearson \
        --whatToPlot scatterplot \
        --plotFile {output.scatterplot} \
        --outFileCorMatrix {params.scatterplot_table} \
        2> {log.LOG_local_log};
        """

#################################################################################
### deepTools: plot coverage information
#################################################################################

rule plot_coverage:
    input:
        bam = expand("{output_dir}/bam_coverages/{sample}_sorted.bam", output_dir=config["output_dir"], sample=get_all_samples()),
        bai = expand("{output_dir}/bam_coverages/{sample}_sorted.bam.bai", output_dir=config["output_dir"], sample=get_all_samples())
    output:
        coverage_plot = "{output_dir}/coverage_plot.png"
    params:
        LOG_cluster_log = "{output_dir}/cluster_log/plot_coverage.log",
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = "{output_dir}/local_log/plot_coverage.log",
    resources:
        threads = 8,
        mem = 20000
    benchmark:
        "{output_dir}/cluster_log/plot_coverage_benchmark.log"
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        plotCoverage \
        --bamfiles {input.bam} \
        --plotFile {output.coverage_plot} \
        --numberOfProcessors {resources.threads} \
        2> {log.LOG_local_log};
        """

#################################################################################
### convert fasta genome to .2bit format
#################################################################################

rule fa_to_2bit:
    input:
        genome = config["genome"]
    output:
        genome_2bit = "{output_dir}/genome.2bit"
    params:
        LOG_cluster_log = "{output_dir}/cluster_log/fa_to_2bit.log",
        queue = "6hours",
        time = "1:00:00"
    log:
        LOG_local_log = "{output_dir}/local_log/fa_to_2bit.log",
    resources:
        threads = 8,
        mem = 20000
    benchmark:
        "{output_dir}/cluster_log/fa_to_2bit_benchmark.log"
    conda:
        "env_yaml/ucsc.yaml"
    shell:
        """
        faToTwoBit {input.genome} {output.genome_2bit} \
        2> {log.LOG_local_log};
        """

#################################################################################
### deepTools: GC bias plots
#################################################################################

rule GC_plot:
    input:
        bam = "{output_dir}/bam_coverages/{sample}_sorted.bam",
        bai = expand("{output_dir}/bam_coverages/{sample}_sorted.bam.bai", output_dir=config["output_dir"], sample=get_all_samples()),
        genome_2bit = "{output_dir}/genome.2bit"
    output:
        gc_plot = "{output_dir}/GC_bias/{sample}_gc.png"
    params:
        gc_tsv = "{output_dir}/GC_bias/{sample}_gc.tsv",
        LOG_cluster_log = "{output_dir}/cluster_log/GC_plot_{sample}.log",
        queue = "6hours",
        time = "6:00:00"
    log:
        LOG_local_log = "{output_dir}/local_log/GC_plot_{sample}.log",
    resources:
        threads = 8,
        mem = 20000
    benchmark:
        "{output_dir}/cluster_log/GC_plot_{sample}_benchmark.log"
    conda:
        "env_yaml/deeptools.yaml"
    shell:
        """
        computeGCBias \
        --bamfile {input.bam} \
        --effectiveGenomeSize 2913022398 \
        --genome {input.genome_2bit} \
        --numberOfProcessors {resources.threads} \
        --GCbiasFrequenciesFile {params.gc_tsv} \
        --biasPlot {output.gc_plot} \
        --plotFileFormat png \
        2> {log.LOG_local_log};
        """




# generate effective genome sizes