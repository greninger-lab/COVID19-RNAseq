#Fig 5 DE by sex as a function of covid status for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu

setwd("~/Desktop")

library(tidyverse)
library(forestmangr)

metadata_file <- file.choose()   
metadata <- read_csv(metadata_file)
all_names <- metadata$alt_name
metadata <- filter(metadata, sex != "U")
metadata <- column_to_rownames(metadata, "alt_name") #name of column containing sample names
samples <- rownames(metadata) 

counts_file <- file.choose() 
counts <- read_csv(counts_file)
counts <- round_df(counts, digits = 0, rf = "round")
counts <- column_to_rownames(counts, "gene") #name of column containing gene names
colnames(counts) <- all_names
counts <- counts[, samples]


#normalization
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design= ~ covid_status + sex + covid_status:sex ) #last is for DEG, earlier are as counfounders (conf1 + conf2 + forDEG)
#This should find differences that result from the interaction between sex and covid status
keep <- (rowSums(counts(dds)) >= 431) #pre-filter to remove any genes without avg 1 count per sample 
#this makes downstream processing much less computationally expensive
dds <- dds[keep,]
dds <- DESeq(dds, parallel=TRUE, quiet = FALSE) 
res <- results(dds)
sum(res$pvalue <0.05, na.rm = TRUE)
sum(res$padj <0.1, na.rm = TRUE)
resSig <- res[which(res$padj < 0.1 ), ] #this is somewhat arbitrary, just to cut down on number of lines. 
resSig <- as.data.frame(resSig)
resSig <- rownames_to_column(resSig)
resSig <- dplyr::arrange(resSig, padj)
norm_counts <- as.data.frame(counts(dds, normalized = TRUE)) #add replaced = TRUE to include replacements for outliers
sig_genes <- resSig$rowname
norm_sig <- norm_counts[sig_genes,]

#assign chromosomes to sig_genes
library("biomaRt")
listMarts()
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
BM <- getBM(attributes=c("external_gene_name", "chromosome_name"),
      filters = "external_gene_name", 
      values = sig_genes, 
      mart = ensembl)
BM2 <- dplyr::filter(BM, !grepl("CHR", chromosome_name))
colnames(BM2)[1] <- "rowname"
resSig <- inner_join(resSig, BM2)

#csv output
write.csv(results(dds), "~/Desktop/MF_results_pre-filtered.csv")
write.csv(resSig, "~/Desktop/MF_significant_genes.csv")
write.csv(as.data.frame(counts(dds, normalized = TRUE)), "~/Desktop/MF_norm_counts_pre-filtered.csv")
write.csv(norm_sig, "~/Desktop/MF_sig_gene_counts.csv") 

