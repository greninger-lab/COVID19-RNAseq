#Fig 5 DE by age as a function of covid status for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu

setwd("~/Desktop")

library(tidyverse)
library(forestmangr)

metadata_file <- file.choose()   
metadata <- read_csv(metadata_file)
all_names <- metadata$alt_name
metadata <- filter(metadata, age != "NA")
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
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design= ~ sixty_plus + covid_status + sixty_plus:covid_status) #last is for DEG, earlier are as counfounders (conf1 + conf2 + forDEG)
#This should test for differences that result from the interaction between age and covid status
dds$sixty_plus <- factor(dds$sixty_plus, levels = c("FALSE","TRUE")) #first position is reference; if not explicit, will be alphabetical
keep <- (rowSums(counts(dds)) >= 467) #pre-filter to remove any genes without avg 1 count per sample 
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

#csv output
write.csv(results(dds), "~/Desktop/age_covid-int_results_pre-filtered.csv")
write.csv(resSig, "~/Desktop/age_covid-int_significant_genes.csv")
write.csv(as.data.frame(counts(dds, normalized = TRUE)), "~/Desktop/age_covid-int_norm_counts_pre-filtered.csv")
write.csv(norm_sig, "~/Desktop/age_covid-int_sig_gene_counts.csv") 

