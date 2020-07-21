#Fig 4 longitudinal analysis for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu


library(tidyverse)
library(forestmangr)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(viridis)
library(cowplot)

setwd("~/Desktop")

metadata_file <- file.choose()   
metadata <- read_csv(metadata_file)
all_samples <- metadata$sample
metadata <- column_to_rownames(metadata, "sample") 
metadata <- metadata[c("Pt_1_coll_1","Pt_1_coll_2", "Pt_2_coll_1", "Pt_2_coll_2", "Pt_7_coll_1", "Pt_7_coll_2"),] #exclude any columns of data that didn't pass QC, either by name or index
sample_names <- rownames(metadata)

counts_file <- file.choose() 
counts <- read_csv(counts_file)
colnames(counts) <- c("gene", all_samples) #rename samples
counts <- counts[,c("gene", sample_names)] #filter samples
counts <- round_df(counts, digits = 0, rf = "round")
counts <- column_to_rownames(counts, "gene") #name of column containing gene names

metadata$collection_number <- as.character(metadata$collection_number) 
metadata$pt <- as.character(metadata$pt) #because patient shouldn't be continuous variable


#Differential expression
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design= ~ pt + collection_number) #last is for DEG, earlier are as counfounders (conf1 + conf2 + forDEG)
keep <- (rowSums(counts(dds)) >= 6) & (rowSums(q) > 2) #pre-filter to remove any genes without avg 1 count per sample AND more than 2 samples with non-zero counts
#this is REALLY aggressive pre-filtering for a small sample set
dds <- dds[keep,]
dds <- DESeq(dds, parallel=FALSE, quiet = FALSE) 
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
write.csv(results(dds), "~/Desktop/results_pre-filtered.csv")
write.csv(resSig, "~/Desktop/significant_genes.csv")
write.csv(as.data.frame(counts(dds, normalized = TRUE)), "~/Desktop/norm_counts_pre-filtered.csv")
write.csv(norm_sig, "~/Desktop/sig_gene_counts.csv")

results <- as.data.frame(res)
results <- rownames_to_column(results)

#add entrez ids, most analyses below rely on these
hs <- org.Hs.eg.db
entrez <- AnnotationDbi::select(hs, 
                                keys = results$rowname,
                                columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
entrez <- distinct(entrez, SYMBOL, .keep_all = TRUE) #removes duplicate entries if there was multiple mapping
results$entrez <- entrez$ENTREZID
results_withentrez <- filter(results, entrez != "NA") #remove entries where there was not mapping to entrez id. this should mostly affect poorly characterized genes 

#Make geneList and gene objects used throughout. Genes are those that are both padj <0.1 and l2fc >1
geneList <- results_withentrez[,3] #L2FC
names(geneList) = as.character(results_withentrez[,8])
sig <- filter
sig <- filter(results_withentrez, padj < 0.1)
sig <- filter(sig, abs(log2FoldChange) > 1)
sig <- arrange(sig, desc(abs(log2FoldChange)))
gene <- sig[,8]


#4B GO enrichment
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",        #BP, CC, or MF
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)
go_df <- as.data.frame(ego@result)
#remove duplicate categories, keeping one with pvalue
go_df <- distinct(go_df, geneID, .keep_all = TRUE)
go_df_20 <- go_df[c(1,4,6,7,13:28),]
des <- go_df_20$Description
des <- rev(des)
go_df_20$Description <- factor(go_df_20$Description, levels = go_df_20$Description[order(-(go_df_20$pvalue))])

png(filename="~/Desktop/4B_long_GO_plot.png", width=5, height=4, units="in", res=300) #fix this, padjust wrong
ggplot(go_df_20, aes(Description, Count, fill=p.adjust)) + 
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = str_wrap(des, width = 40)) +
  scale_fill_viridis() +
  xlab("GO Term") +
  ylab("Number Enriched") +
  theme(
    legend.position = c(0.82,0.3),
    text = element_text(size = 10),
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10))
dev.off()

#####4C genes

plot <- column_to_rownames(results_withentrez, "rowname") 

hum <- go_df_20[5, "geneID"]
hum <- as.vector(str_split_fixed(hum, "/", 11))
hum <- str_sort(hum)
H <- plot[hum,]
H <- rownames_to_column(H)

humoral <- ggplot(H, aes(x=rowname, y = log2FoldChange)) +
  geom_col(color="black", fill= "dark grey") +
  geom_errorbar( aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width = .25) +
  geom_hline(yintercept=0) +
  ylab("Log 2 Fold Change") +
  ylim(-7,9) +
  ggtitle("Humoral Immune Response") +
  theme(
    title = element_text(size = 9),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size=7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7),
    legend.position = "none")

wound <- go_df_20[15, "geneID"]
wound <- as.vector(str_split_fixed(wound, "/", 17))
wound <- str_sort(wound)
W <- plot[wound,]
W <- rownames_to_column(W)

wound <- ggplot(W, aes(x=rowname, y = log2FoldChange)) +
  geom_col(color="black", fill= "dark grey") +
  geom_errorbar( aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width = .25) +
  geom_hline(yintercept=0) +
  ylab("Log 2 Fold Change") +
  ylim(-7,9) +
  ggtitle("Wound Healing") +
  theme(
    title = element_text(size = 9),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size=7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7),
    legend.position = "none")

covid <- c("CXCL10", "IFIT3", "IFIT2", "CXCL9", "IFIT1", "RSAD2", "OASL", "GBP1", "IFI44L", "OAS3", "CXCL11", "HERC5", "DDX58", "IFI44", "MX1", "OAS1", "UBE2L6", "OAS2", "HERC6")
results <- column_to_rownames(results, "X")
C <- plot[covid,]
C <- rownames_to_column(C)
C$rowname <- factor(C$rowname, levels = C$rowname[order(C$padj)])
C$threshold <- C$padj < 0.1

consensus <- ggplot(C, aes(x=rowname, y = log2FoldChange, fill=threshold)) +
  geom_col(color="black") +
  geom_errorbar( aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width = .25) +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("white", "dark grey")) +
  ylab("Log 2 Fold Change") +
  ylim(-8,3) +
  ggtitle("SARS-CoV-2 Consensus Genes") +
  theme(
    title = element_text(size = 9),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size=7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7),
    legend.position = "none")

rib <- go_df_20[2, "geneID"]
rib <- as.vector(str_split_fixed(rib, "/", 24))
R <- plot[rib,]
R <- rownames_to_column(R)

ribosome <- ggplot(R, aes(x=rowname, y = log2FoldChange)) +
  geom_col(color="black", fill= "dark grey") +
  geom_errorbar( aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width = .25) +
  geom_hline(yintercept=0) +
  ylab("Log 2 Fold Change") +
  ylim(0,8) +
  ggtitle("Ribosomal Proteins") +
  theme(
    title = element_text(size = 9),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size=7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7),
    legend.position = "none")

png(filename = "~/Desktop/4C_long_genes.png", width = 9, height = 5, units = "in", res=300)
cowplot::plot_grid(humoral, wound, consensus, ribosome, ncol=2)
dev.off()



