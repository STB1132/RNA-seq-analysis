# RNA-seq-analysis

This repository contains a **Nextflow pipeline** to analyze RNA-Seq data from SRA samples.  
The pipeline performs the following steps:


1. Perform quality control with FastQC
2. Align reads to the reference genome with HISAT2
3. Count reads per gene with featureCounts
4. Merge counts and perform differential expression analysis with DESeq2
5. Extract expression of FOXP3 and identify most differentially expressed genes

But first, to get the samples, you can use [fasterqdump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) or [kingfisher](https://wwood.github.io/kingfisher-download/). Here I recomend using kingfisher since its easier to use, faster and gives less errors. Why using this isntead of an usual wget or download? Files can get corrupted due to many different reasons (fluctuations of wifi conecction, errors in the dowload) and the dowloaded file can be corrupted. So, as a good standard, please use any of the resourse availiable like this one. To dowload the fastq files just get [their docker image](https://hub.docker.com/r/wwood/kingfisher/tags). Then run:

```bash
docker run --rm \
  -v /path/to/your/foolder:/data \
  -w /data \
  wwood/kingfisher:0.4.1 \
  get -r SRR5223500 \ // Or whatever SRR id you want
  -m aws-http prefetch \
  --output-format-possibilities fastq.gz \
  -t 4

```


---

## Requirements

- [Docker](https://www.docker.com/products/docker-desktop/)
- [Nextflow](https://www.nextflow.io/)
- Mac M1 / M2 (ARM64) or Linux

---

## Setup

### 1. Build Docker Image

For Mac M1/M2, use the following command to build the Docker image with `amd64` emulation:

```bash
docker buildx build --platform linux/amd64 -t rna_seq_pipeline .
RNA-seq data quality filter, standardisation and analysis using scalability and reproducibility standards. 
```
Alternatevely, you can run build the docker image as usual if you dont have the ARM64 aquitecture

 ```bash
 docker  build  -t rna_seq_pipeline .
```
