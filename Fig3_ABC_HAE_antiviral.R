#Fig 3 intrinsic antiviral consensus genes for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu

setwd("~/Desktop")

library(tidyverse)
library(VennDiagram)
library(SuperExactTest)
library(pheatmap)

#read in gene lists generated earlier

HAE_100 <- read.csv("~/Desktop/Figures/HAE_top100DE.csv", stringsAsFactors = FALSE) 
posneg_83 <- read.csv("~/Desktop/Figures/posneg_all83DE.csv", stringsAsFactors = FALSE) 
hilo_all <- read.csv("~/Desktop/Figures/hilo_all_sig.csv", stringsAsFactors = FALSE)

posneg_top83DE <- posneg_83$gene
hilo_allDE <- hilo_all$gene
HAE_top100DE <- moscona_100$gene

venn_all <- list(posneg_top83DE, hilo_allDE, HAE_top100DE)

venn.diagram(venn_all, 
             category.names = c("Positive vs Negative" , "High vs Low Viral Load" , "Infected vs Uninfected HAE"),
             filename = 'venn_allDE.png',
             output=FALSE,
             
             # Output features
             imagetype="png", height = 900, width = 900, resolution = 300, compression = "lzw",
             # Circles
             lwd = 2, lty = 'blank', fill = c("red", "green", "blue"),
             # Numbers
             cex = .6, fontface = "bold", fontfamily = "sans",
             # Set names
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer", cat.pos = c(-27, 27, 180), cat.dist = c(0.055, 0.055, 0.055), cat.fontfamily = "sans", rotation = 1)

list <- list(posneg_top83DE, hilo_allDE, HAE_top100DE)
results <- supertest(list, n=35272)

#heatmap
l2fc <- read_csv("~/Desktop/consensus_genes_L2FC.csv") #csv containing l2fc from consensus genes for HAE, pos/neg, hi/low
l2fc <- column_to_rownames(l2fc, "X1")

png(filename="~/Desktop/3X_l2fc_consensus.png", width=2, height=3, units="in", res=300)
print(pheatmap(l2fc, show_colnames = FALSE, fontsize_row = 6, cluster_cols=FALSE, cluster_rows = FALSE)) #creates dendrogram for samples and genes. if cluster_cols=FALSE it does not cluster
dev.off()


###3B
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)
library(viridis)

covid <- c("CXCL10", "IFIT3", "IFIT2", "CXCL9", "IFIT1", "RSAD2", "OASL", "GBP1", "IFI44L", "OAS3", "CXCL11", "HERC5", "DDX58", "IFI44", "MX1", "OAS1", "UBE2L6", "OAS2", "HERC6")

hs <- org.Hs.eg.db
entrez <- AnnotationDbi::select(hs, 
                                 keys = covid,
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "SYMBOL")
entrez <- distinct(entrez, SYMBOL, .keep_all = TRUE) #removes duplicate entries if there was multiple mapping
covid <- as.data.frame(covid)
covid$entrez <- entrez$ENTREZID

edgn <- enrichDGN(covid$entrez, readable = TRUE) #disease networks 
res <- as.data.frame(edgn@result)
res20 <- res[1:20,]
res20["umls:C0375023", "Description"] <- "Respiratory Syncytial Virus"
des <- res20$Description
des <- rev(des)
res20$Description <- factor(res20$Description, levels = res20$Description[order(-(res20$p.adjust))])

png(filename="~/Desktop/3B_DGN_plot.png", width=5, height=4, units="in", res=300)
ggplot(res20, aes(Description, Count, fill=p.adjust)) + 
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = str_wrap(des, width = 50)) +
  scale_fill_viridis() +
  ylab("Number Enriched") +
  xlab("Disease") +
  theme(
    legend.position = c(0.8,0.3),
    text = element_text(size = 10),
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10))
dev.off()


 #3C Network
PN_results <- read.csv("~/Desktop/Figures/posneg_results_pre-filtered.csv")
colnames(covid)[1] <- "gene"
colnames(PN_results)[1] <- "gene"
covid2 <- inner_join(covid, PN_results)

covidList <- covid2[,4] #L2FC
names(covidList) = as.character(covid2[,2])

png(filename="~/Desktop/3C_DGN_networked_genes.png", width=5, height=4, units="in", res=300)
cnetplot(edgn, categorySize="geneNum", foldChange=covidList)
dev.off()

