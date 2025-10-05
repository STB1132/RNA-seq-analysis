#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
count_file <- args[1]
sample_file <- args[2]

library(DESeq2)
library(readr)
library(dplyr)

counts <- read_delim(count_file, "\t")
samples <- read_delim(sample_file, "\t")

count_matrix <- counts %>% select(-Geneid)
rownames(count_matrix) <- counts$Geneid

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = samples %>% select(sample, condition) %>% column_to_rownames("sample"),
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

# Save results
write.csv(as.data.frame(res), "results/deseq2/deseq2_results.csv")

# FOXP3 expression
foxp3 <- plotCounts(dds, gene="FOXP3", intgroup="condition", returnData=TRUE)
write.csv(foxp3, "results/deseq2/FOXP3_expression.csv", row.names=FALSE)

