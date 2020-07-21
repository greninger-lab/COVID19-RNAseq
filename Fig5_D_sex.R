##Fig 5D sex differences for Lieberman et al, PLoS Biology 2020
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


#5D

res$`Transcript Type` <- c("Non-immune", "B cell", "Non-immune", "Non-immune", "Non-immune", "Negative Immune Regulator", "Non-immune", "NK Activating Receptor", "Non-immune", "Negative Immune Regulator", "Non-immune", "Non-immune", "Non-immune", "Negative Immune Regulator", "B cell", "Non-immune", "Negative Immune Regulator", "Non-immune", "Non-immune")

png(filename="~/Desktop/5Dsexdiff_barplot.png", width=8, height=5, units="in", res=300)
ggplot(res, aes(x= reorder(SYMBOL, -log2FoldChange), y = log2FoldChange, fill=`Transcript Type`)) +
  geom_col(color="black") +
  coord_flip() +
  geom_errorbar( aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width = .25) +
  geom_hline(yintercept=0) +
  scale_fill_manual(values = c("#ffcc00", "#ff0000", "#33cc00", "dark grey")) +
  xlab("gene") +
  ylab("Log 2 Fold Change") +
  ylim(-8,6) +
  theme(
    legend.position = c(0.8,0.8),
    axis.text.y = element_text(size = 10, face = "italic"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 8)
    )
dev.off()

