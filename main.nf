#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samples = "samples.tsv"
params.genome  = "GRCh38" 
params.gtf     = "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"
params.fasta   = "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

process DOWNLOAD_SRA {
    tag "$sample_id"
    container "quay.io/biocontainers/sra-tools:3.0.6--pl5321h9f5acd7_0"
    publishDir "results/fastq", mode: 'copy'

    input:
    tuple val(sample_id), val(sra), val(condition)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz")

    script:
    """
    fasterq-dump $sra --split-files --gzip -e 4
    mv ${sra}_1.fastq.gz ${sample_id}_1.fastq.gz
    mv ${sra}_2.fastq.gz ${sample_id}_2.fastq.gz
    """
}

process FASTQC {
    tag "$sample_id"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    publishDir "results/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(reads1), path(reads2)

    output:
    tuple val(sample_id), val(condition), path(reads1), path(reads2), path("fastqc_reports")

    script:
    """
    mkdir fastqc_reports
    fastqc -t 4 -o fastqc_reports $reads1 $reads2
    """
}

process ALIGN_HISAT2 {
    tag "$sample_id"
    container "quay.io/biocontainers/hisat2:2.2.1--h1b792b2_5"
    publishDir "results/bam", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(reads1), path(reads2)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}.bam")

    script:
    """
    hisat2-build -p 4 ${params.fasta} index
    hisat2 -p 4 -x index -1 $reads1 -2 $reads2 | samtools sort -@4 -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    """
}

process FEATURECOUNTS {
    tag "$sample_id"
    container "quay.io/biocontainers/subread:2.0.6--h7132678_1"
    publishDir "results/counts", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(bam)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}.counts.txt")

    script:
    """
    featureCounts -T 4 -a ${params.gtf} -o ${sample_id}.counts.txt $bam
    """
}

process MERGE_COUNTS {
    publishDir "results/merged", mode: 'copy'
    container "rocker/tidyverse:4.3.1"

    input:
    tuple val(sample_id), val(condition), path(count_files)

    output:
    path("merged_counts.txt")

    script:
    """
    Rscript -e '
    library(dplyr)
    files <- list.files("results/counts", pattern="counts.txt", full.names=TRUE)
    counts <- Reduce(function(x, y) merge(x, y, by="Geneid"), lapply(files, function(f) read.delim(f, comment.char="#")))
    write.table(counts, "merged_counts.txt", sep="\\t", quote=FALSE, row.names=FALSE)
    '
    """
}

process DESEQ2_ANALYSIS {
    publishDir "results/deseq2", mode: 'copy'
    container "rocker/tidyverse:4.3.1"

    input:
    path("merged_counts.txt")

    script:
    """
    Rscript bin/deseq2_analysis.R merged_counts.txt ${params.samples}
    """
}

workflow {
    Channel
        .fromPath(params.samples)
        .splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.sample, row.sra, row.condition) }
        |>
        DOWNLOAD_SRA
        |>
        FASTQC
        |>
        ALIGN_HISAT2
        |>
        FEATURECOUNTS
        |>
        MERGE_COUNTS
        |>
        DESEQ2_ANALYSIS
}

