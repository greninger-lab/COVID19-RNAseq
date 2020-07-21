##Fig S2C of adjusted bacterial RPM by viral load for Lieberman et al, PLoS Biology 2020
#naliebe@uw.edu

library(tidyverse)

setwd("~/Desktop/")

meta <- read_csv("~/Desktop/Figures/metadata_all_for_figs.csv") 
bac <- read_csv("~/Desktop/Figures/adj_rpm.csv") #contains the virus-subtracted bacterial RPM

plot <- inner_join(bac, meta, by = "acc_num")
plot <- filter(plot, Viral_Load != "Unknown")

plot <- plot[c(13:49, 51:463),]
plot$Viral_Load <- factor(plot$Viral_Load, levels = c("Negative", "Low", "Mid", "High"))


png(filename="~/Desktop/S2c_bacterialRPM.png", width=2, height=3, units ="in", res = 300)
ggplot(plot, aes( x=Viral_Load, y=bacterial_rpm, fill=Viral_Load)) +
  geom_boxplot(outlier.size=1) + 
  geom_jitter(color="black", size=0.4, alpha = 0.5, width = 0.25) +
  scale_y_log10() +
  scale_x_discrete(labels = str_wrap(c("Negative", "Low (N1 ct > 24)", "Mid (N1 ct 24-19)", "High (N1 ct < 19)"), width = 9)) +
  scale_fill_manual(values = c("#ff0000", "#ffcc00", "#ff9900", "#33cc00")) +
  ylab("Bacterial RPM") +
  xlab("Viral Load") +
  theme(
    legend.position = "none",
    title = element_text(size = 7),
    text = element_text(size = 5),
    axis.title.x = element_text(size = 6), 
    axis.title.y = element_text(size = 6))
dev.off()

#Stats
#between groups
neg <- filter(plot, plot$Viral_Load == "Negative")
low <- filter(plot, plot$Viral_Load == "Low")
mid <- filter(plot, plot$Viral_Load == "Mid")
high <- filter(plot, plot$Viral_Load == "High")

#normal?
shapiro.test(neg$bacterial_rpm) #p-value = 9.723e-10 NOT NORMAL
shapiro.test(low$bacterial_rpm) #p-value < 2.2e-16
shapiro.test(mid$bacterial_rpm) #p-value < 2.2e-16
shapiro.test(high$bacterial_rpm) #p-value < 2.2e-16

wilcox.test(neg$bacterial_rpm, low$bacterial_rpm) #7.141e-07
wilcox.test(mid$bacterial_rpm, low$bacterial_rpm) #0.8838
wilcox.test(high$bacterial_rpm, low$bacterial_rpm) #0.006688
wilcox.test(neg$bacterial_rpm, mid$bacterial_rpm) #1.172e-07
wilcox.test(neg$bacterial_rpm, high$bacterial_rpm) #0.0002822
wilcox.test(mid$bacterial_rpm, high$bacterial_rpm) #0.001413

#across groups

kruskal.test(bacterial_rpm ~ Viral_Load, data = plot) #0.001413

