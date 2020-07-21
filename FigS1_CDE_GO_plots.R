##Fig Supp 1C,D,E GO Analysis for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(viridis)

setwd("~/Desktop/")

results <- read.csv("~/Desktop/Figures/posneg_results_pre-filtered.csv", stringsAsFactors = FALSE)

#add entrez ids
hs <- org.Hs.eg.db
entrez <- AnnotationDbi::select(hs, 
                                keys = results$X,
                                columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
entrez <- distinct(entrez, SYMBOL, .keep_all = TRUE) #removes duplicate entries if there was multiple mapping
results$entrez <- entrez$ENTREZID
results <- filter(results, entrez != "NA") #remove entries where there was not mapping to entrez id. this should mostly affect poorly characterized genes 

geneList <- results[,3] #L2FC
names(geneList) = as.character(results[,8])
sig <- filter(results, padj < 0.1)
sig <- filter(sig, abs(log2FoldChange) > 1)
sig <- arrange(sig, desc(abs(log2FoldChange)))
gene <- sig[,8]

#GO enrichment
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",        #BP, CC, or MF
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)
go_df <- as.data.frame(ego@result)
#remove duplicate categories, keeping one with lowest pvalue
go_df <- distinct(go_df, geneID, .keep_all = TRUE)
go_df_20 <- go_df[1:20,]

des <- go_df_20$Description
des <- rev(des)
go_df_20$Description <- factor(go_df_20$Description, levels = go_df_20$Description[order(-(go_df_20$pvalue))])
#go_df_20$Count <- factor(go_df_20$Count, levels = go_df_20$Count)

#S1C

png(filename="~/Desktop/S1C_GO_plot.png", width=5, height=6, units="in", res=300)
ggplot(go_df_20, aes(Description, Count, fill=p.adjust)) + 
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = str_wrap(des, width = 30)) +
  scale_fill_viridis() +
  ylab("Number Enriched") +
  xlab("GO Term") +
  theme(
    legend.position = "right",
    text = element_text(size = 10),
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10))
dev.off()

#S1D

res_table <- filter(results, padj != "NA")
abs <- abs(res_table$log2FoldChange)
res_table$abs <- abs
res_table_ordered <- arrange(res_table, desc(abs))
res_table_ordered <- column_to_rownames(res_table_ordered, "X")

imm <- go_df_20[2, "geneID"]
imm <- as.vector(str_split_fixed(imm, "/", 20))
B <- res_table_ordered[imm,]

png(filename="~/Desktop/S1D_viralimmune_plot.png", width=5, height=3, units="in", res=300)
ggplot(B, aes(x=rownames(B), y = log2FoldChange)) +
  geom_col(color="black", fill="dark grey") +
  geom_errorbar( aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width = .25) +
  geom_hline(yintercept=0) +
  ylim(-0.5,5.5) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, face = "italic")) + 
  xlab("Gene") +
  ylab("Log2 Fold Change")
dev.off()


#S1E

rib <- go_df_20[1, "geneID"]
rib <- as.vector(str_split_fixed(rib, "/", 16))
C <- res_table_ordered[rib,]

png(filename="~/Desktop/S1E_RP_plot.png", width=5, height=3, units="in", res=300)
ggplot(C, aes(x=rownames(C), y = log2FoldChange)) +
  geom_col(color="black", fill="dark grey") +
  geom_errorbar( aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width = .25) +
  geom_hline(yintercept=0) +
  ylim(-3, 0.5) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, face = "italic")) + 
  xlab("Gene") +
  ylab("Log2 Fold Change")
dev.off()
