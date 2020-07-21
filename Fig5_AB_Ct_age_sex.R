##Fig 5A-B age and sex by ct for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu

library(tidyverse)
library(cowplot)

setwd("~/Desktop/")

alldata <- read.csv("~/Desktop/Figures/meta_plus_normcounts.csv", stringsAsFactors = FALSE)
bac <- read.csv("~/Desktop/Figures/metadata_for_boxplots.csv")
alldata <- full_join(bac, alldata, by ="acc_num")
alldata <- filter(alldata, Viral_Load != "Unknown")
alldata <- filter(alldata, age_bin != "unk")

adj <- (alldata[,22:16055]) + 1 #If values are zero, can't be plotted properly
adj <- cbind(alldata[1:20], adj)

pos_adj <- filter(adj, covid_status == "pos")
pos_adj <- filter(pos_adj, sex != "U")
neg_adj <- filter(adj, covid_status == "neg")

#colors <- c("#9900cc", "#990066", "#990000")

#boxplot 5A Ct by Age
png(filename="~/Desktop/Ct_by_age_bin.png", width=3, height=2, units ="in", res = 300)
ggplot(pos_adj, aes( x=age_bin, y=covid_N1_ct, fill=age_bin)) +
  geom_boxplot(outlier.size=1) + 
  #geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  geom_point(size = 0.3, position = position_jitterdodge(jitter.width = 0.4)) +
  #scale_y_log10() +
  scale_fill_manual(values = c("#9900cc", "#990066", "#990000")) +
  ylab("N1 Ct") +
  xlab("Age") +
  #ggtitle("CCL2") +
  theme(
    legend.position = "none",
    title = element_text(size = 8),
    text = element_text(size = 6),
    axis.title.x = element_text(size = 7), 
    axis.title.y = element_text(size = 7))
dev.off()

#boxplot 5B Ct by sex
png(filename="~/Desktop/Ct_by_sex.png", width=2, height=2, units ="in", res = 300)
ggplot(pos_adj, aes( x=sex, y=covid_N1_ct, fill=sex)) +
  geom_boxplot(outlier.size=1) + 
  #geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  geom_point(size = 0.3, position = position_jitterdodge(jitter.width = 0.4)) +
  #scale_y_log10() +
  scale_fill_manual(values = c("#66cc66","#33ccff")) +
  ylab("N1 Ct") +
  xlab("Sex") +
  #ggtitle("CCL2") +
  theme(
    legend.position = "none",
    title = element_text(size = 8),
    text = element_text(size = 6),
    axis.title.x = element_text(size = 7), 
    axis.title.y = element_text(size = 7))
dev.off()

#Stats
pos <- filter(alldata, covid_status == "pos")
neg <- filter(alldata, covid_status == "neg")

young <- filter(pos, age_bin == "0-29")
mid <- filter(pos, age_bin == "30-59")
old <- filter(pos, age_bin == "60-100")

M <- filter(pos, sex == "M")
F <- filter(pos, sex == "F")

shapiro.test(F$covid_N1_ct) #mid: p=0.03 so need to use non parametric. M and F normal

#across groups
kruskal.test(covid_N1_ct ~ age_bin, data = pos) #aov for normally distributed data
#p for pos$age_bin= 0.66

#between groups
t.test(M$covid_N1_ct, F$covid_N1_ct) 
#p=0.883

