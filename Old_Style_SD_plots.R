#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")


sedef=readbed("../Assembly_analysis/SEDEF/{V}.SDs.bed", "T2T CHM13")
sedef=sedef[fracMatch >= 0.95 & alnB >= 10000]
inter= sedef[chr2!=chr]
intra = sedef[chr2==chr & start < start2]
dim(sedef)
dim(inter)
dim(intra)
N=100
if(F){
intra= intra[sample(.N,N)]
inter= inter[sample(.N,N)]
}

h = 9/1.5
w = 16/1.5
#intra$color = GRAY
#intra.new = overlap_either(intra, NEW)
#intra$color[intra.new] = NEWCOLOR
#intra$color[intra.new & intra$fracMatch >= 0.95 & intra$alnB >= 10000] = "#ff0000"
intra$color = NEWCOLOR
intra = intra[order(-color)]
print(intra %>% group_by(color) %>% summarise(length(chr)))


kpinit <- function(){
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, chromosomes = CHRS[1:length(CHRS)-1])
  return(kp)
}
kpintra <-  function(kp){
  kpPlotLinks(kp, 
              data=toGRanges(data.frame(chr=intra$chr,start=intra$start,intra$end) ), data2=toGRanges(data.frame(chr=intra$chr2,start=intra$start2,intra$end2) ), 
              col="darkblue", border = NA)
}
kpinter <- function(kp){
  kpPlotLinks(kp, 
              data=toGRanges(data.frame(chr=inter$chr,start=inter$start,inter$end) ), data2=toGRanges(data.frame(chr=inter$chr2,start=inter$start2,inter$end2)), 
              col=transparent(NEWCOLOR, amount=0), border = NA)
}

Intra=as.ggplot(expression(
  kp<-kpinit(),
  kpintra(kp)
)); 
ggsave("old_intra.pdf", plot = Intra, height = h, width = w)

Inter=as.ggplot(expression(
  kp<-kpinit(),
  kpinter(kp)
)); 
ggsave("old_inter.pdf", plot = Inter, height = h, width = w)

Both=as.ggplot(expression(
  kp<-kpinit(),
  kpinter(kp),
  kpintra(kp),
  kpAddCytobands(kp)
))
ggsave("old_both.pdf", plot = Both, height = h, width = w)
  
Bothr=as.ggplot(expression(
  kp<-kpinit(),
  kpintra(kp),
  kpinter(kp),
  kpAddCytobands(kp)
))
ggsave("old_bothr.pdf", plot = Bothr, height = h, width = w)

