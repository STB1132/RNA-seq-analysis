#!/usr/bin/env nextflow


nextflow.enable.dsl=2


params.samples = "/Users/stb1132/Desktop/RNA-seq-analysis/samples.tsv"
params.genome  = "GRCh38"
params.gtf     = "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"
params.fasta   = "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" 
params.fasta_gz_local = "${projectDir}/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
params.fasta_local = "${projectDir}/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.index_dir = "s3://genome-idx/hisat/grch38_genome.tar.gz"


process PREPARE_HISAT2_INDEX {
    container "ubuntu:22.04"

    input:
    path index_tar

    output:
    path "hisat2_index", emit: index_dir

    script:
    """
    mkdir hisat2_index
    tar -xzf ${index_tar} -C hisat2_index --strip-components=1
    """
}

process DOWNLOAD_SRA {
    tag "$sample_id"

    container "quay.io/biocontainers/sra-tools:3.1.0--h9f5acd7_0"

    input:
    tuple val(sample_id), val(condition), val(sample_path)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz")

    script:
    """
    prefetch ${sample_id} -p -v
    fasterq-dump --split-files --threads 4  --outdir . "${sample_id}"
    gzip -f "${sample_id}_1.fastq" "${sample_id}_2.fastq"

    """
}

process FASTQC {
    tag "$sample_id"

    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    publishDir "results/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(read1), path(read2)

    output:
    tuple val(sample_id), val(condition), path(read1), path(read2)

    script:
    """
    fastqc -t 4 $read1 $read2
    """
}

process ALIGN_HISAT2 {
    tag "$sample_id"

    container "nanozoo/hisat2:2.2.1.commit7e01700--5e923e8"

    publishDir "results/bam", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(read1), path(read2)
    path index_dir

    output:
    tuple val(sample_id), val(condition), path("${sample_id}.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.bam.bai"), emit: bai

    script:
    """
    hisat2 \
        -p 4 \
        -x ${index_dir}/genome \
        -1 $read1 \
        -2 $read2 \
    | samtools sort -@ 4 -o ${sample_id}.bam -

    samtools index ${sample_id}.bam
    """
}


workflow {

    // Prepare HISAT2 index from S3 tarball
    index_ch = Channel.fromPath(params.index_dir)
    hisat2_index = PREPARE_HISAT2_INDEX(index_ch)

    reads = Channel
        .fromPath(params.samples)
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            tuple(
                row.sample,
                row.condition,
                row.path
            )
        }

    download_out = DOWNLOAD_SRA(reads)

    qc_out = FASTQC(download_out)

    bam_out = ALIGN_HISAT2(qc_out, hisat2_index.index_dir)

    bam_out.bam.view()
}