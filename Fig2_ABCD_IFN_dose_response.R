#Fig 2 interferon dose response for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu

library(tidyverse)
library(forestmangr)
library(DESeq2)
library(cowplot)
library(ggrepel)


alldata <- read.csv("~/Desktop/Figures/meta_plus_normcounts.csv", stringsAsFactors = FALSE)
alldata <- filter(alldata, Viral_Load != "Unknown")
alldata$Viral_Load <- factor(alldata$Viral_Load, levels = c("Negative", "Low", "Mid", "High"))

#2A: IFN gene boxplots by dose
#adding 1 to all counts to allow plotting by log scale without visual artifacts. 
#ANOVA on data with zero

adj <- (alldata[,21:16054]) + 1
adj <- cbind(alldata[1:20], adj)

colors <- c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")

ACE2 <- ggplot(adj, aes( x=Viral_Load, y=ACE2, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("ACE2 ****") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_text(size = 6, face = "plain"))

TMPRSS2 <- ggplot(adj, aes( x=Viral_Load, y=TMPRSS2, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("TMPRSS2") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_blank())

CXCL9 <- ggplot(adj, aes( x=Viral_Load, y=CXCL9, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("CXCL9 ****") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_blank())

OASL <- ggplot(adj, aes( x=Viral_Load, y=OASL, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  #ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("OASL ****") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_blank())

MX1 <- ggplot(adj, aes( x=Viral_Load, y=MX1, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("MX1 ****") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_blank())

CD274 <- ggplot(adj, aes( x=Viral_Load, y=CD274, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  #ylab("Normalized Counts") +
  xlab("Viral Load") +
  ylab("Normalized Counts") +
  ggtitle("CD274/PD-L1 ****") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_text(size = 6, face = "plain"))

USP18 <- ggplot(adj, aes( x=Viral_Load, y=USP18, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  #ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("USP18 ****") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_blank())

CCL2 <- ggplot(adj, aes( x=Viral_Load, y=CCL2, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  #ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("CCL2 ****") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_blank())

RPL4 <- ggplot(adj, aes( x=Viral_Load, y=RPL4, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("RPL4") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_blank())

RPS6 <- ggplot(adj, aes( x=Viral_Load, y=RPS6, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("RPS6") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_blank())

png(filename="~/Desktop/2A_IFN.png", width=8, height=3, units="in", res=300)
cowplot::plot_grid(ACE2, TMPRSS2, CXCL9, OASL, MX1, CD274, USP18, CCL2, RPL4, RPS6,  ncol = 5, nrow = 2) 
dev.off()

#mann whitney test since none are normal
neg <- filter(alldata, alldata$Viral_Load == "Negative")
low <- filter(alldata, alldata$Viral_Load == "Low")
mid <- filter(alldata, alldata$Viral_Load == "Mid")
high <- filter(alldata, alldata$Viral_Load == "High")


wilcox.test(high$ACE2, low$ACE2) #t.test replaces wilcox.text for gaussian data

#2B data prep

metadata_file <- file.choose()   
metadata <- read_csv(metadata_file)
hilo_metadata <- filter(metadata, covid_status == "pos")
hilo_metadata <- hilo_metadata[!is.na(hilo_metadata$hi_v_lo),]
metadata <- column_to_rownames(metadata, "alt_name") #name of column containing sample names
hilo_metadata <- column_to_rownames(hilo_metadata, "alt_name")
hilo_samples <- rownames(hilo_metadata)

counts_file <- file.choose() 
counts <- read_csv(counts_file)
counts <- round_df(counts, digits = 0, rf = "round")
counts <- column_to_rownames(counts, "gene") #name of column containing gene names
colnames(counts) <- rownames(metadata)
hilo_counts <- counts[, hilo_samples]

#normalization
dds <- DESeqDataSetFromMatrix(countData = hilo_counts, colData = hilo_metadata, design= ~ seq_run + hi_v_lo) #last is for DEG, earlier are as counfounders (conf1 + conf2 + forDEG)
dds$hi_v_lo <- factor(dds$hi_v_lo, levels = c("hi","lo")) #first position is reference; this refers to hi ct samples (= low viral load)
keep <- (rowSums(counts(dds)) >= 207)  #pre-filter to remove any genes without avg 1 count per sample 
dds <- dds[keep,]
dds <- DESeq(dds, parallel=FALSE, quiet = FALSE) 
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
res <- results(dds)
sum(res$pvalue <0.05, na.rm = TRUE)
sum(res$padj <0.01, na.rm = TRUE)
resSig <- res[which(res$padj < 0.1 ), ]
resSig <- as.data.frame(resSig)
resSig <- rownames_to_column(resSig)
resSig <- dplyr::arrange(resSig, padj)
resFold <- dplyr::arrange(resSig, desc(abs(log2FoldChange)))

sig_genes <- resSig$rowname
norm_sig <- norm_counts[sig_genes,]

#csv output
write.csv(results(dds), "~/Desktop/hilo_results_pre-filtered.csv")
write.csv(resSig, "~/Desktop/hilo_significant_genes.csv")
write.csv(as.data.frame(counts(dds, normalized = TRUE)), "~/Desktop/hilo_norm_counts_pre-filtered.csv")
write.csv(norm_sig, "~/Desktop/hilo_sig_gene_counts.csv") 


###2B lo vs hi viral load volcano
hilo <- as.data.frame(res)
res_table <- filter(hilo, padj != "NA")
sig <- filter(res_table, padj < 0.1)
abs <- abs(res_table$log2FoldChange)
res_table$abs <- abs
threshold <- (res_table$abs > 1.5) & (res_table$padj < 0.05)
res_table$threshold <- threshold

sig <- arrange(sig, by = desc(log2FoldChange))
sig$genelabels <- ""
sig$genelabels[1:15] <- sig$X[1:15] #15 highest l2fc up with 
sig$genelabels[186:200] <- sig$X[186:200] #15 lowest l2fc down

res_table_ordered <- arrange(res_table, desc(abs))
res_table_ordered <- full_join(res_table_ordered, sig, by = "X")
res_table_ordered <- column_to_rownames(res_table_ordered, "X")

png(filename="~/Desktop/2B_volcano_15up_15down.png", width=5, height=5, units="in", res=300)
ggplot(res_table_ordered, aes(x=log2FoldChange.x, y=-log10(padj.x), color=threshold, label = genelabels)) +
  geom_point() +
  geom_text_repel(color = "black", size = 2.5, box.padding = 0.3, segment.size = 0.3) +
  geom_vline(xintercept = -1.5, linetype = "dotted", color = "grey30", size = 0.25) +
  geom_vline(xintercept = 1.5, linetype = "dotted", color = "grey30", size = 0.25) +
  geom_hline(yintercept = 1.301, linetype = "dotted", color = "grey30", size = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  xlab("Log2 Fold Change") + 
  ylab("-Log10 Adjusted p-Value") +
  xlim(-5, 6)  +
  ylim(-1,21) +
  theme(
    legend.position = "none",
    text = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )
dev.off()

#2C CIBERSORT

meta <- read.csv("~/Desktop/Figures/metadata_all_for_figs.csv", stringsAsFactors = FALSE)
meta <- filter(meta, covid_status == "pos")
meta <- filter(meta, Viral_Load != "Unknown")
meta <- filter(meta, Viral_Load != "Mid")
meta <- filter(meta, pass_manual_QC == "TRUE")

col <- str_c(meta$Viral_Load, meta$alt_name, sep="_")

mat <- read_csv("~/Desktop/Figures/CIBERSORT_byHiLo.csv")
mat$Mixture <- col

mat <- mat[, 1:23]
mat <- arrange(mat, desc(col))
mat <- column_to_rownames(mat, "Mixture")
t_mat <- as.data.frame(t(mat))
t_mat <- rownames_to_column(t_mat)
t_mat$Low <- rowMeans(t_mat[,2:100])
t_mat$High <- rowMeans(t_mat[,101:208])

av <- t_mat[,c(1,209,210)]
av$rowname <- factor(av$rowname, levels = c("Neutrophils**", "Eosinophils", "Mast cells activated**", "Mast cells resting", "Dendritic cells activated**", "Dendritic cells resting", "Macrophages M2*", "Macrophages M1****", "Macrophages M0", "Monocytes", "NK cells activated*", "NK cells resting**", "T cells gamma delta*", "T cells regulatory (Tregs)", "T cells follicular helper", "T cells CD4 memory activated", "T cells CD4 memory resting", "T cells CD4 naive**", "T cells CD8", "Plasma cells", "B cells memory", "B cells naive**"))
colnames(av)[1] <- "Cell Type"
av_melt <- melt(av, id.vars="Cell Type")

colors <- c("black", "grey30", "grey50", "grey70", "#990066", "#990099", "#990000", "#cc0000", "#cc0066", "#cc3399", "#3399ff", "#33ccff", "#339900", "#33cc00", "#33cc66", "#669900", "#66cc00", "#66cc99","#66ff00", "#cc9900", "#cccc00", "#ffcc00")

png("~/Desktop/2Ca_cibersort_AVct.png", width=2, height=5, unit = "in", res=300)
ggplot(av_melt, aes(x=variable, y=value, fill=`Cell Type`)) +
  geom_col() + 
  scale_fill_manual(values = colors) +
  xlab("Viral Load") +
  ylab("Proportion") +
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8))
dev.off()

png("~/Desktop/2Cb_legend.png", width=5.1, height=5, unit = "in", res=300)
ggplot(av_melt, aes(x=variable, y=value, fill=`Cell Type`)) +
  geom_col() + 
  scale_fill_manual(values = colors) +
  xlab("Viral Load") +
  ylab("Proportion") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 12))
dev.off()

#joined image and legend in illustrator

#2D cell types

IGHA1 <- ggplot(adj, aes( x=Viral_Load, y=IGHA1, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("IGHA1 *") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_text(size = 6, face = "plain"))

CD22 <- ggplot(adj, aes( x=Viral_Load, y=CD22, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("CD22 *") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_blank())

CXCL8 <- ggplot(adj, aes( x=Viral_Load, y=CXCL8, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("CXCL8 **") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_text(size = 6, face = "plain"))

S100A9 <- ggplot(adj, aes( x=Viral_Load, y=S100A9, fill=Viral_Load)) +
  geom_violin(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#33cc00", "#ffcc00", "#ff9900", "#ff0000")) +
  ylab("Normalized Counts") +
  xlab("Viral Load") +
  ggtitle("S100A9 *") +
  theme(
    legend.position = "none",
    title = element_text(size = 7, face = "italic"),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6, face = "plain"), 
    axis.title.y = element_blank())

png(filename="~/Desktop/2D_B_neut.png", width=3.2, height=3, units="in", res=300)
cowplot::plot_grid(IGHA1, CD22, CXCL8, S100A9,  ncol = 2, nrow = 2) 
dev.off()


wilcox.test(high$S100A9, low$S100A9) 


