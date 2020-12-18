#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
library(RColorBrewer)
#library(dplyr)
library(grid)
#library(gridBase)
library(gridExtra)
library(data.table)
library(gtable)
#source("http://bioconductor.org/biocLite.R")
#biocLite("karyoploteR")
#BiocManager::install("karyoploteR")
library(karyoploteR)
library(GenomicRanges)
library(cowplot)
library(glue)
library(viridis)
library(colourvalues)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#suppressPackageStartupMessages(library("argparse"))
#library(argparse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
CHRS= c(paste0("chr",seq(1,22)),"chrX")


in_rgn <-function(df, rgn){
  #df = tbl_nof
  #rgn = "chr6:28,510,120-33,480,577"
  rgn = gsub(",", "", rgn)
  x = unlist(strsplit(rgn, "[:-]"))
  chr = x[1]
  start = as.integer(x[2])
  end = as.integer(x[3])
  cond = (chr == df$reference_name) & (start <= df$reference_end & df$reference_start <= end)
  #print(sum(cond))
  return(cond)
}

tbl = fread("../Masker/split_alignments/{V}.split.nosd.bed"); tbl$RGN = "Non SD"
tbl2 = fread("../Masker/split_alignments/{V}.split.sd.bed"); tbl2$RGN = "SD"
tbl3 = fread("../Masker/split_alignments/{V}.split.sdflank.bed"); tbl3$RGN = "SD_Flanks"
tbl_nof = rbind(tbl, tbl2, tbl3)
colnames(tbl_nof)[1]=sub("#", "", colnames(tbl_nof)[1]) 


# fix the queyrt coords 
tbl_nof=separate(tbl_nof,query_name, "-|:", into=c("query_name","q_add","q_trash"))
tbl_nof$q_add = as.numeric(tbl_nof$q_add) -1 
tbl_nof$query_start = tbl_nof$q_add + tbl_nof$query_start
tbl_nof$query_end = tbl_nof$q_add + tbl_nof$query_end


tbl_nof$chr = tbl_nof$query_name
tbl_nof$chr = factor(tbl_nof$chr, 
                                levels =  c(CHRS, unique(tbl_nof$chr[which(!tbl_nof$chr %in% CHRS)]) ) , ordered = TRUE)
tbl_nof=tbl_nof[order(chr)]


tbl_nof$Divergence = 100-tbl_nof$perID_by_events
tbl_nof$aln_frac =  (tbl_nof$query_end - tbl_nof$query_start) / tbl_nof$query_length

chrX = tbl_nof[tbl_nof$reference_name == "chrX"]
chrX$RGN="chrX"
tbl_nof = rbind(tbl_nof, chrX)


mhc_rgn = "chr6:28,510,120-33,480,577"
mhc = tbl_nof[ in_rgn(tbl_nof, mhc_rgn) ]
mhc$RGN = "MHC"
tbl_nof = rbind(tbl_nof, mhc)


colors = c(`SD`="#af0404", `MHC`="#3282b8", `chrX`="#96bb7c", `SD_Flanks`="#ede682", `Non SD`="#000000")


#tbl = tbl_nof[tbl_nof$perID_by_events > quantile(tbl_nof$perID_by_events, 0.025)]
tbl=tbl_nof
# add colors 
ncol = 10
ii <- cut(tbl$perID_by_events, breaks = seq(min(tbl$perID_by_events), max(tbl$perID_by_events), len = ncol), include.lowest = TRUE)
#breaks=unique(c(quantile(tbl$perID_by_events, probs = seq(0, 1, by = 1/(ncol-1) ) ))); 
#ncol = length(breaks)+1
#ii <- cut(tbl$perID_by_events, breaks=breaks, include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
#tbl$color = colorRampPalette(inferno(ncol-1, direction=-1, end=0.85) )(ncol-1)[ii]
tbl$color = rev(colorRampPalette(c("black", "red"), alpha=F, bias=.1)(ncol-1))[ii]


p = ggplot(data = tbl, aes(Divergence)) +
  geom_histogram(aes(fill=color), alpha=0.75, bins = 50)+
  facet_grid(RGN ~ ., scales = "free_y") + 
  scale_fill_identity() +      
  #scale_y_log10() + annotation_logticks(sides="l")+
  theme_cowplot(); p

p2 = ggplot(data = tbl_nof, aes(Divergence+.01, color = RGN)) + stat_ecdf(size=1.5, alpha=0.75) + 
  scale_x_log10() + annotation_logticks(sides="b")+
  scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
  xlab("Percent divergence") +
  ylab("Cumulative fraction of 5kbp windows") +
  theme_classic()+
  theme(legend.position = "top") + guides(fill=guide_legend(ncol=length(colors))) ; p2

p3 = ggplot(data=tbl) + geom_hex(aes(x=aln_frac, y=perID_by_events), bins=50) + 
  scale_fill_viridis(name = "count", breaks = 10^(0:10),  trans = "log")+
  facet_grid(~RGN) + 
  theme_cowplot();
ggsave("short_alns.pdf", plot = p3, height = 4, width = 16)


grid = plot_grid(p2) + ggtitle("Divergence of 5kbp windows of CHM13 aligned to GRCh38"); grid
ggsave("SD_divergence.pdf", plot=grid,  height = 4.5, width = 8)
ggsave("SD_divergence.png", plot=grid,  height = 4.5, width = 8, dpi=600)
grid







if(T){
pdf(file = "DivergedIdeogram.pdf", height = 9, width = 9)
SD = tbl[tbl$RGN=="SD" & tbl$perID_by_events < 99.5]
genome = toGRanges(data.frame(tbl_nof %>% group_by(chr) %>% summarise(start=0, end=max(query_end))))
kp = plotKaryotype(genome = genome)
kpPoints(kp, chr=SD$chr, x=SD$query_start, y=SD$Divergence,col=SD$color, ymax = 20)
#kpAxis(kp, ymax = 15)
#kpPlotDensity(kp, range,data.panel = 2, window.size=2e6, ymax =1000)
kpAddMainTitle(kp,"Diverged CHM13 SDs (<99.5)")
kpAddBaseNumbers(kp)
dev.off()

}
if(T){
  pdf(file = "DivergedIdeogram_hg38.pdf", height = 9, width = 16)
  SD = tbl[tbl$RGN=="SD" & tbl$perID_by_events < 99.5]
  #genome = toGRanges(data.frame(tbl_nof %>% group_by(chr) %>% summarise(start=0, end=max(query_end))))
  kp = plotKaryotype(genome = "hg38")
  kpPoints(kp, chr=SD$reference_name, x=SD$reference_start, y=SD$Divergence, col=SD$color, ymax = 20)
  #kpAxis(kp, ymax = 15)
  #kpPlotDensity(kp, range,data.panel = 2, window.size=2e6, ymax =1000)
  kpAddMainTitle(kp,"Diverged CHM13 SDs (<99.5)")
  kpAddBaseNumbers(kp)
  dev.off()
}

if(F){
  genes =fread("../Masker/split_alignments/{V}.split.sd.genes.bed", col.names = c("chr", "start","end","name","perid"));

kp = plotKaryotype(genome="hg38", plot.type = 2)
SD = tbl[tbl$RGN=="SD" & tbl$perID_by_events < 99.5]
sd_genes = genes[genes$perid < 99.5]; sd_genes$Divergence = 100-sd_genes$perid
setorder(SD, Divergence)
range = makeGRangesFromDataFrame(SD[,c("reference_name","reference_start","reference_end")], 
                                 seqnames.field = "reference_name", start.field = "reference_start", end.field = "reference_end")
#genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
#sd_genes = subsetByOverlaps(gene, range, ignore.strand=T)

#kpPlotRegions(kp, range, col = SD$color, avoid.overlapping = F, border = NA)
kpPlotMarkers(kp, chr = sd_genes$chr, x=sd_genes$start, y=0.5, labels = sd_genes$name, cex=0.25, line.color="gray")

kpPoints(kp, chr=SD$reference_name, x=SD$reference_start, y=SD$Divergence,col=SD$color, ymax = 10)
kpAxis(kp, ymax = 5)
kpPlotDensity(kp, range,data.panel = 2, window.size=2e6, ymax =1000)
kpAddMainTitle(kp,"Diverged CHM13 SDs (<99.5)")
kpAddBaseNumbers(kp)
dev.off()
write.table(sd_genes, row.names = F, 
            file="sd_diverged_genes.txt", sep = "\t", quote = F)
}

sum( (tbl$query_end - tbl$query_start)[tbl$RGN=="SD"])/10^6



separate(tbl_nof,query_name, "-|:", into=c("query_name","q_add","q_trash"))
