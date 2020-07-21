##Fig S2B of adjusted bacterial RPM by viral load for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu

library(tidyverse)
library(pheatmap)

setwd("~/Desktop/")

metadata <- read_csv("~/Desktop/metagenomic_heatmap.csv") #contains RPM values for viral and bacterial species of interest
metadata <- column_to_rownames(metadata, "alt_name")
metadata2 <- arrange(metadata, by = desc(covid_N1_ct))

#make annotation dataframe, the part we are interested in
df <- metadata2[, c(18, 3, 23, 30:35)] #indices for columns about to be named, below
colnames(df) <- c("Viral Load", "Covid Status", "Other Viruses", "S. aureus", "H. influenzae", "P. aeruginosa", "S. pneumoniae", "S. pseudopneumoniae", "M. catarrhalis")

###only reading the count data in because pheatmap doesn't seem to be able to just generate the annotation data, cropping that portion out instead.
PN_results <- read.csv("~/Desktop/posneg_results_pre-filtered.csv", stringsAsFactors = FALSE)
PN_norm <- read.csv("~/Desktop/posneg_norm_counts_pre-filtered.csv", stringsAsFactors = FALSE)
PN_norm <- column_to_rownames(PN_norm, "X")
#for heatmap, but won't use the heatmap portion
by_ct <- rownames(metadata2)
PN_norm <- PN_norm[, by_ct]
PN_results <- arrange(PN_results, padj)
PN_sig <- filter(PN_results, padj <0.1)
PN_sig <- arrange(PN_sig, by = desc(abs(log2FoldChange)))
PN_sig <- filter(PN_sig, abs(log2FoldChange) > 1) %>% arrange(padj) #original, for fig 1
PN_50sig <- PN_sig$X[1:50]
PN_de <- (PN_norm[PN_50sig,]) + 1
fc <- PN_de / rowMeans(PN_de)
l2fc <- log2(fc)
l2fc <- l2fc[, by_ct] #again, these values don't matter here


rownames(df) <- colnames(l2fc) #don't matter but need to match

my_colour = list(
  `Covid Status` = c(neg = "#33cc00", pos = "#0066ff"),
  `Viral Load` = c(Negative = "#33cc00", Low = "#ffcc00", Mid = "#ff9900", High = "#ff0000", Unknown = "white"),
  `Other Viruses` = c(None = "black", Enterovirus = "#33ffff", Enterovirus_and_RSV = "#33cccc", InfluenzaB = "#ccff00", Metapneumovirus = "#cc99ff", HKU1 = "#ff00ff", NL63 = "#ff0099"),
  `S. aureus` = c(D = "black", C = "grey50", B = "grey80", A = "yellow"),
  `H. influenzae` = c(D = "black", C = "grey50", B = "grey80", A = "yellow"),
  `P. aeruginosa` = c(D = "black", C = "grey50", B = "grey80", A = "yellow"),
  `S. pneumoniae` = c(D = "black", C = "grey50", B = "grey80", A = "yellow"),
  `S. pseudopneumoniae` = c(D = "black", C = "grey50", B = "grey80", A = "yellow"),
  `M. catarrhalis` = c(D = "black", C = "grey50", B = "grey80", A = "yellow"))
  

png(filename="~/Desktop/S2B_heatmap.png", width=15, height=12, units="in", res=300)
print(pheatmap(l2fc, annotation_col=df, annotation_colors = my_colour, fontsize = 10, fontsize_row = 6, treeheight_col = 10, treeheight_row = 10, show_colnames = FALSE, cluster_cols=FALSE))
dev.off()
