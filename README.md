# RNA-seq-analysis
# RNA-Seq Analysis Pipeline

This repository contains a **Nextflow pipeline** to analyze RNA-Seq data from SRA samples.  
The pipeline performs the following steps:

1. Download SRA samples
2. Perform quality control with FastQC
3. Align reads to the reference genome with HISAT2
4. Count reads per gene with featureCounts
5. Merge counts and perform differential expression analysis with DESeq2
6. Extract expression of FOXP3 and identify most differentially expressed genes

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
