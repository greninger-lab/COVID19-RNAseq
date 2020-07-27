##Fig S3B sex differences for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(viridis)

mf <- read.csv("~/Desktop/Figures/MF_significant_genes.csv") #from DE of sex differences
counts <- read.csv("~/Desktop/Figures/posneg_norm_counts_pre-filtered.csv")
counts <- column_to_rownames(counts, "X")

hs <- org.Hs.eg.db
entrez <- AnnotationDbi::select(hs, 
                                keys = rownames(counts),
                                columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
entrez <- distinct(entrez, SYMBOL, .keep_all = TRUE) #removes duplicate entries if there was multiple mapping
entrez <- filter(entrez, ENTREZID != "NA")

colnames(mf)[2] <- "SYMBOL"
res <- inner_join(mf, entrez, by = "SYMBOL")

geneList <- entrez$ENTREZID
gene <- res$ENTREZID

ego <- enrichGO(gene          = gene,
                universe      = geneList,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",        #BP, CC, or MF
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)
go_df <- as.data.frame(ego@result)
go_df <- distinct(go_df, geneID, .keep_all = TRUE)

go_df_10 <- go_df[1:10,]
des <- go_df_10$Description
des <- rev(des)
go_df_10$Description <- factor(go_df_10$Description, levels = go_df_10$Description[order(-(go_df_10$p.adjust))])

png(filename="~/Desktop/S3B_sexdiff_GO_plot.png", width=5, height=4, units="in", res=300)
ggplot(go_df_10, aes(Description, Count, fill=p.adjust)) + 
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = str_wrap(des, width = 30)) +
  scale_fill_viridis() +
  ylab("Number Enriched") +
  xlab("GO Term") +
  theme(
    legend.position = c(0.85,0.25),
    text = element_text(size = 10),
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10))
dev.off()
