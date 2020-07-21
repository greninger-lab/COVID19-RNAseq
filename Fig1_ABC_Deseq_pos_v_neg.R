#DEseq and figure 1 for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu

setwd("~/Desktop")

library(tidyverse)
library(forestmangr)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(viridis)


metadata_file <- file.choose()   
metadata <- read_csv(metadata_file)
metadata <- column_to_rownames(metadata, "alt_name") 
all_names <- rownames(metadata)

counts_file <- file.choose() 
counts <- read_csv(counts_file)
counts <- round_df(counts, digits = 0, rf = "round")
counts <- column_to_rownames(counts, "gene") 
colnames(counts) <- all_names

#normalization
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design= ~ seq_run + covid_status) #last is for DEG, earlier are as counfounders (conf1 + conf2 + forDEG)
keep <- (rowSums(counts(dds)) >= 484) #pre-filter to remove any genes without avg 1 count per sample 
dds <- dds[keep,]
dds <- DESeq(dds, parallel=TRUE, quiet = FALSE) 
res <- results(dds)
sum(res$pvalue <0.05, na.rm = TRUE)
sum(res$padj <0.1, na.rm = TRUE)
resSig <- res[which(res$padj < 0.1 ), ]  
resSig <- as.data.frame(resSig)
resSig <- rownames_to_column(resSig)
resSig <- dplyr::arrange(resSig, padj)
norm_counts <- as.data.frame(counts(dds, normalized = TRUE)) 
sig_genes <- resSig$rowname
norm_sig <- norm_counts[sig_genes,]

#csv output
write.csv(results(dds), "~/Desktop/posneg_results_pre-filtered.csv")
write.csv(resSig, "~/Desktop/posneg_significant_genes.csv")
write.csv(as.data.frame(counts(dds, normalized = TRUE)), "~/Desktop/posneg_norm_counts_pre-filtered.csv")
write.csv(norm_sig, "~/Desktop/posneg_sig_gene_counts.csv") 

#Figure 1A: heatmap of top 50 most significant genes. 
df <- metadata[, c("Viral_Load", "covid_status")]  #add columns as desired to add to metadata. 
colnames(df) <- c("Viral Load", "Covid Status")
de <- (norm_sig[1:50,] + 1)
fc <- de / rowMeans(de)
l2fc <- log2(fc)

png(filename="~/Desktop/1A_heatmap.png", width=90, height=20, units="in", res=300)
print(pheatmap(l2fc, annotation_col=df, cluster_cols=TRUE)) #creates dendrogram for samples and genes. if cluster_cols=FALSE it does not cluster
dev.off()

##Figure 1B Volcano plot of top 15 highest abs(l2FC)
PN_results <- as.data.frame(res)
res_table <- filter(PN_results, padj != "NA")
abs <- abs(res_table$log2FoldChange)
res_table$abs <- abs
threshold <- (res_table$abs > 1.5) & (res_table$padj < 0.05)
res_table$threshold <- threshold

PN_sig <- arrange(PN_sig, by = desc(log2FoldChange))
PN_sig$genelabels <- ""
PN_sig$genelabels[1:15] <- PN_sig$X[1:15] #15 highest l2fc up
PN_sig$genelabels[69:83] <- PN_sig$X[69:83] #15 lowest l2fc down

res_table_ordered <- arrange(res_table, desc(abs))
res_table_ordered <- full_join(res_table_ordered, PN_sig, by = "X")
res_table_ordered <- column_to_rownames(res_table_ordered, "X")


png(filename="~/Desktop/1B_volcano.png", width=5, height=3, units="in", res=300)
ggplot(res_table_ordered, aes(x=log2FoldChange.x, y=-log10(padj.x), color=threshold, label = genelabels)) +
  geom_point() +
  geom_text_repel(color = "black", size = 2.5, box.padding = 0.5, segment.size = 0.3) +
  geom_vline(xintercept = -1.5, linetype = "dotted", color = "grey30", size = 0.5) +
  geom_vline(xintercept = 1.5, linetype = "dotted", color = "grey30", size = 0.5) +
  geom_hline(yintercept = 1.301, linetype = "dotted", color = "grey30", size = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  xlab("Log2 Fold Change") + 
  ylab("-Log10 Adjusted p-Value") +
  xlim(-5, 5)  +
  theme(
    legend.position = "none",
    text = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )
dev.off()

##Figure 1c GSEA - input made with GSEA software run from desktop
gsea <- read.csv("~/Desktop/Figures/GSEA_NES.csv")
gsea <- arrange(gsea, desc(NES))

gsea$NAME <- factor(gsea$NAME, levels = gsea$NAME[order(gsea$NES)])

png("~/Desktop/1C_gsea.png", width = 5, height = 3, units = "in", res=300)
ggplot(gsea, aes(x=NAME, y=NES, fill=FDR)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis() +
  labs(x = "Gene Set", y = "Normalized Enrichment Score") +
  theme(
    text = element_text(size = 8),
    plot.title = element_blank(),
    legend.title = element_text(size = 10),
    legend.position = c(0.85,0.35),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )
dev.off()
