#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if(! require("tidyverse")) install.packages("tidyverse")
if(! require("ggnewscale")) install.packages("ggnewscale")
if(! require("ggrepel")) install.packages("ggrepel")
if(! require("data.table")) install.packages("data.table")
if(! require("glue")) install.packages("glue")
if(! require("RColorBrewer")) install.packages("RColorBrewer")
if(! require("scales")) install.packages("scales")
if(! require("cowplot")) install.packages("cowplot")


adf = SEDEF_V1.1[acro & acro2]


a = ggplot(data=adf, aes(y=matchB, x=fracMatch, color=intra))+
  geom_point(alpha=0.25, size=0.5)+
  geom_density_2d()+
  scale_y_log10(label=comma)+ annotation_logticks()+
  scale_color_manual(values=c(NEWCOLOR, OLDCOLOR)) +
  ylab("length") + xlab("Identity of the SD alignment") +
  theme_cowplot() + ggtitle("Acrocentric SDs distributed by identity and length")

b = ggplot(data=adf, aes(x=fracMatch, fill=intra))+
  geom_histogram(breaks=seq(0.9, 1, 0.001))+
  scale_fill_continuous()+
  ylab("") + xlab("Identity of the SD alignment") +
  scale_fill_manual(values=c(NEWCOLOR, OLDCOLOR)) +
  theme_cowplot() + ggtitle("# of SD alignments for a given identity")

c = ggplot(data=adf, aes(x=fracMatch , weight=matchB, fill=intra))+
  geom_histogram(breaks=seq(0.9, 1, 0.001))+
  #scale_y_log10(label=comma)+ annotation_logticks()+
  scale_fill_manual(values=c(NEWCOLOR, OLDCOLOR)) +
  ylab("") + xlab("Identity of the SD alignment") +
  theme_cowplot() + ggtitle("# of bp in SD alignments for a given identity") 

p = plot_grid(a,b,c, ncol=1, align = "v", rel_heights = c(2,1,1))
ggsave(glue("{SUPP}/acrop_SDs_length_and_identity.pdf"), plot=p, height = 12, width = 12)
p