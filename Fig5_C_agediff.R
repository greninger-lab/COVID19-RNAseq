##Fig 5C age differences for Lieberman et al, PLoS Biology 2020
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

#Fig 5C
CXCL11 <- ggplot(age, aes( x=age_bin, y=CXCL11, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  #geom_point(size = 0.25, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("CXCL11 ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_text(size = 6, face = "plain"),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

CXCL9 <- ggplot(age, aes( x=age_bin, y=CXCL9, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  #geom_point(size = 0.25, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("CXCL9 ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

CXCL10 <- ggplot(age, aes( x=age_bin, y=CXCL10, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  #geom_point(size = 0.25, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("CXCL10 ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

CXCR3 <- ggplot(age, aes( x=age_bin, y=CXCR3, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  #geom_point(size = 0.25, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("CXCR3 *") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_text(size = 6, face = "plain"),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

GZMB <- ggplot(age, aes( x=age_bin, y=GZMB, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  #geom_point(size = 0.25, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("GZMB ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

CD8A <- ggplot(age, aes( x=age_bin, y=CD8A, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  #geom_point(size = 0.25, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("CD8A ****") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 6),
    title = element_text(size = 7, face = "italic"),
    axis.title.x = element_text(size = 6, face = "plain"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5)
  )

#only for legend, cropping graph itself
ACTB <- ggplot(age, aes( x=age_bin, y=ACTB, fill=`Covid Status`)) +
  geom_boxplot(outlier.size=0.5) + 
  #geom_point(size = 0.25, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_y_log10() +
  scale_fill_manual(values = c("#33cc00", "#0066ff")) +
  xlab("Age") +
  ylab("Normalized Counts") +
  ggtitle("ACTB") +
  theme(
    legend.position = "left",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)
  )

png(filename = "~/Desktop/5C_age.png", width = 6, height = 3, units = "in", res = 300)
cowplot::plot_grid(CXCL11, CXCL9, CXCL10, ACTB, CXCR3, GZMB, CD8A, ncol=4, nrow=2)
dev.off()

