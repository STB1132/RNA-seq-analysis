FROM rocker/r-ver:4.3.1

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    apt-utils \
    sra-toolkit \
    fastqc \
    hisat2 \
    samtools \
    subread \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN Rscript -e "install.packages('tidyverse', repos='https://cloud.r-project.org')" && \
    Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    Rscript -e "BiocManager::install(c('DESeq2'))"

WORKDIR /workspace










