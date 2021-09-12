#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")


sedef=readbed("../Assembly_analysis/SEDEF/{V}.SDs.bed", "T2T CHM13")
inter= sedef[chr2!=chr]
intra = sedef[chr==chr2 & start < start2]
#intra=sample_n(intra,100)
#inter=sample_n(inter,100)

chrs = CHRS[1:length(CHRS)-1]
intra$color = GRAY
intra.new <- overlap_either(intra, NEW)
intra$color[intra.new] = NEWCOLOR
intra$color[intra.new & intra$fracMatch >= 0.95 & intra$alnB >= 10000] = "#ff0000"
intra <- intra[order(-color)]
print(intra %>% group_by(color) %>% summarise(length(chr)))

#
# regions for genomic instablity
#
nbp = function(gi){
  s = start(ranges(reduce(gi)))
  e = end(ranges(reduce(gi)))
  
  print()
}

gi.new.df = intra[intra.new ]%>% filter( fracMatch >= 0.95 & alnB >=10000 & abs(start-start2) < 10e6 & abs(start-start2) > 50000 ) %>% dplyr::select(chr, start, end2)
write.table(gi.new.df, file = "GenomicInstabilityRegions_new.bed", sep="\t", quote = F, row.names = F, col.names = F)
gi.new = toGRanges(gi.new.df )
reduce(gi.new)
sum(width(reduce(gi.new)))/1e6

gi.df=intra%>% filter( fracMatch >= 0.95 & alnB >=10000 & abs(start-start2) < 10e6 & abs(start-start2) > 50000 ) %>% dplyr::select(chr, start, end2)
gi = toGRanges(gi.df)
reduce(gi)
sum(width(reduce(gi)))/1e6

gi.old <- setdiff(reduce(gi), reduce(gi.new))
#
#
#




Intra=as.ggplot(expression(
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, chromosomes = NOM),
  kpPlotRegions(kp, data=reduce(gi.new), col=transparent("orange", 0.25), data.panel = "ideogram", r0=0.0, r1=1, border = NA),
  kpPlotRegions(kp, data=gi.old, col=transparent("orange", 0.25), data.panel = "ideogram", r0=0, r1=1, border = NA),
  
  kpPlotLinks(kp, 
              data=toGRanges(data.frame(chr=intra$chr,start=intra$start,intra$end)), data2=toGRanges(data.frame(chr=intra$chr2,start=intra$start2,intra$end2)), 
              col=intra$color, border = NA)
  #kpAddCytobands(kp)
  #pAddBaseNumbers(kp)
))
ggsave("Intra.pdf", plot=Intra, height=6, width=10)


if(F){
Inter=as.ggplot(expression(
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, chromosomes = chrs),
  kpPlotLinks(kp, 
              data=toGRanges(data.frame(chr=inter$chr,start=inter$start,inter$end) ), data2=toGRanges(data.frame(chr=inter$chr2,start=inter$start2,inter$end2)), 
              col=transparent("darkred", amount=0), border = NA),
  kpAddCytobands(kp),
  kpAddBaseNumbers(kp)
));
}




source("plotutils.R")
hg38 = readbed("data/misc_files/hg38.chr_only.SDs.bed.gz", "GRCh38", chrfilt=TRUE)
hg37 = readbed("data/misc_files/hg19.no_alt.SDs.bed.gz", "GRCh37")
celera = readbed("data/misc_files/Celera_WGSA.SDs.bed.gz", "Celera WGSA")

multiple_sds = rbind(sedef, hg38,hg37[1])#,celera)
sheall = pairwise_plot(multiple_sds)
sheinter = pairwise_plot(multiple_sds[chr!=chr2])
sheintra = pairwise_plot(multiple_sds[chr==chr2])

# title(main="Intrachromosomal SDs")
C1 = plot_grid(sheall$identity, sheall$length, labels = "C", nrow=2) + 
  ggtitle("All SDs")+ theme(plot.title = element_text(hjust = 0.5)); C1
C2 = plot_grid(sheintra$identity, sheintra$length, nrow=2)+ 
  ggtitle("Intrachromosomal SDs") + theme(plot.title = element_text(hjust = 0.5))
C3 = plot_grid(sheinter$identity, sheinter$length, nrow=2) +
  ggtitle("Interchromosomal SDs")+ theme(plot.title = element_text(hjust = 0.5))
She = plot_grid(C1,C2,C3, ncol = 3)
C1

Cmain = plot_grid(sheall$identity, sheall$length+ theme(legend.position = "none"), labels = "C", nrow=1); Cmain #legend(x = "bottom", fill = transparent(c(NEWCOLOR, OLDCOLOR),0.25), legend = c("SDs more than 0.5% diverged from GRCh38", "Density of all SDs")); Cmain

ideo = plot_grid(Intra, NA,  labels = c("a","b"), rel_widths = c(3,2), ncol=2); 

p = plot_grid(ideo, She, nrow=2)
f1 = plot_grid(ideo, Cmain, nrow=2)
  
ggsave("FigureS1.pdf", plot=p, height = 16, width = 16)
ggsave("figures/Figure1.pdf", plot=f1, height = 8, width = 12)

ggsave("PairWise.pdf", plot=C1, height=8, width=8)
