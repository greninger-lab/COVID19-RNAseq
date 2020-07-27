##Fig S3A age differences for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu

library(tidyverse)
library(cowplot)

setwd("~/Desktop/")

data <- read.csv("~/Desktop/Figures/meta_plus_normcounts.csv", stringsAsFactors = FALSE) # merged metadata with normalized counts
age_stats <- filter(data, age_bin != "unk") #use unadjusted matrix for stats
adj <- age_stats[,21:16054] + 1 #add 1 count to 
age <- cbind(age_stats[1:20], adj)
colnames(age)[4] <- "Covid Status"

age$Viral_Load <- factor(age$Viral_Load, levels = c("Negative", "Low", "Mid", "High"))


#Fig S3A
PCGF6 <- ggplot(age, aes( x=age_bin, y=PCGF6, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("PCGF6") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7 , face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_text(size = 6, face = "plain"),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

ACE2 <- ggplot(age, aes( x=age_bin, y=ACE2, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("ACE2 ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

TMPRSS2 <- ggplot(age, aes( x=age_bin, y=TMPRSS2, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("TMPRSS2") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

DDX58 <- ggplot(age, aes( x=age_bin, y=DDX58, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("DDX58 ***") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

GBP1 <- ggplot(age, aes( x=age_bin, y=GBP1, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("GBP1 ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

HERC5 <- ggplot(age, aes( x=age_bin, y=HERC5, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("HERC5 *") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_text(size = 6, face = "plain"),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

HERC6 <- ggplot(age, aes( x=age_bin, y=HERC6, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("HERC6 ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

IFI44 <- ggplot(age, aes( x=age_bin, y=IFI44, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("IFI44 **") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

IFI44L <- ggplot(age, aes( x=age_bin, y=IFI44L, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("IFI44L ***") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

IFIT1 <- ggplot(age, aes( x=age_bin, y=IFIT1, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("IFIT1 ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

IFIT2 <- ggplot(age, aes( x=age_bin, y=IFIT2, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("IFIT2 **") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_text(size = 6, face = "plain"),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

IFIT3 <- ggplot(age, aes( x=age_bin, y=IFIT3, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("IFIT3 ***") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

MX1 <- ggplot(age, aes( x=age_bin, y=MX1, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("MX1 **") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

OAS1 <- ggplot(age, aes( x=age_bin, y=OAS1, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("OAS1 ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

OAS2 <- ggplot(age, aes( x=age_bin, y=OAS2, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("OAS2 ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )


OAS3 <- ggplot(age, aes( x=age_bin, y=OAS3, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("OAS3 ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_text(size = 6, face = "plain"),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

OASL <- ggplot(age, aes( x=age_bin, y=OASL, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("OASL ***") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

RSAD2 <- ggplot(age, aes( x=age_bin, y=RSAD2, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("RSAD2") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

UBE2L6 <- ggplot(age, aes( x=age_bin, y=UBE2L6, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("UBE2L6 ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

png(filename = "~/Desktop/S3A_age.png", width = 8, height = 6, units = "in", res = 300)
cowplot::plot_grid(PCGF6, ACE2, TMPRSS2, DDX58, GBP1, HERC5, HERC6, IFI44, IFI44L, IFIT1, IFIT2, IFIT3, MX1, OAS1, OAS2, OAS3, OASL, RSAD2, UBE2L6, ACTB, ncol=5)
dev.off()



