#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")

chm13 = readbed("../wssd/kmer/bed/chm13_t2t_wssd.sort.bed","T2T CHM13");chm13
hg38 = readbed("../wssd/kmer/bed/hg38_wssd.sort.bed","GRCh38");hg38

df = bind_rows(chm13,hg38)
df

ggplot(data=sample_n(df,10000), aes(x=y, fill=Assembly)) + geom_histogram(position = "identity", binwidth=1, alpha=0.5) +
  scale_y_continuous(trans="log10") + coord_cartesian(xlim=c(0,10))
  

  
  
  p=ggplot(data=sample_n(df[chr!="chrMT" & V10 < 0.5 ],10000))+
    geom_histogram(aes(x=start/1000000, weight=(end-start), fill=Assembly), binwidth = 10, alpha=0.5)+
    facet_wrap(chr~., ncol=2)+
    geom_segment(aes(x=0,xend=max(chrlen/1000000),y=0,yend=0))+
    scale_x_continuous(labels = comma)+
    xlab("Genomic position (Mbp)")+
    theme_cowplot()+theme(legend.position = "bottom")+
    theme(strip.text = element_text(size = 10, margin = margin()), strip.background = element_blank());p
  