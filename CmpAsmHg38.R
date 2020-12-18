#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)
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
library(ggplotify)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(RColorBrewer)

# check is rstudio and if so set the working direcotry to curdir 
if(rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
GAPS = "../assemblies/ucsc.hg38.no_alts.fasta.gap.bed"
gaps = fread(GAPS, col.names = c("chr", "start", "end"))


#suppressPackageStartupMessages(library("argparse"))
#library(argparse)

read_dupmasker <- function(f, paf=NULL){
  dups = fread(f, col.names = c("chr","start","end","name","score","strand","thickStart","thickEnd","color") )
  dups$hex = sapply(strsplit(dups$color, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255)) 
  
  
  if(!is.null(paf)){
    chrs = unique(paf[,c("chr","len")])
    dups = merge(dups, chrs[,c("chr", "len")], by="chr")
    
    #major_strand = paf %>% group_by(chr, strand) %>% summarise( total = sum(end-start)) %>% 
    #  group_by(chr) %>% summarise(flip = strand[which.max(total)])
    #dups = merge(dups, major_strand[, c("chr", "flip")], by="chr", all.x=T)
  }
  return(dups)
}

read_rm <- function(f){
  rm = fread(f)
  names(rm)[1] = "chr"
  rm = rm[ rm$repeat_class %in% c("Satellite", "Simple_repeat", "Low_complexity") ]
  rm$color = "#cd707b"
  rm$color[rm$repeat_class == "Simple_repeat"] = "#90EE90"
  rm$color[rm$repeat_class == "Low_complexity"] = "#90EE90"
  return(rm)
}

read_fai <-function(f){
  asm = fread(f, col.names = c("chr","end","offset","width","width2" ))
  asm$start=1
  asm = asm[,c("chr", "start", "end")]
  asm$len = asm$end
  return(asm)
}

read_paf <- function(f){
  f = PAF
  paf = fread(f,  
              fill=T,
              col.names=c("chr","len","start","end", "strand", "t_chr","t_len", "t_start", "t_end", "match", "length","qual"),
              drop=13:50)
  # get primary chromosome of alignment 
  prichr <- paf %>% group_by(chr, t_chr) %>% summarise( total = sum(end-start)) %>%
    group_by(chr) %>% summarise(prichr=t_chr[which.max(total)]) 
  paf = merge(paf, prichr, by="chr", all.x=T)
  return(paf)
}

read_diff <- function(f){
  diff = fread(f, fill=TRUE, col.names = c("chr", "type", "start", "end", "lenn", "qlenn", "qdiff"))
  diff = diff[, c("chr", "start", "end", "type",  "lenn", "qlenn", "qdiff") ];  
  diff = diff[order(diff$type),]
  diff$start = diff$start - 1
  
  col_vector = c("#66cdaa", "#ffa500", "#00ff00", "#0000ff", "#1e90ff", "#ff1493")
  names(col_vector) = c("GAP", "DUP", "BRK", "JMP", "INV", "SEQ")
  
  diff$color <- col_vector[match(diff$type, names(col_vector))]; 
  colortbl = unique(diff[, c("type", "color")]); n = nrow(colortbl);
  colortbl=colortbl[order(colortbl$type),]
  colortbl = unique(diff[, c("type", "color")]); n = nrow(colortbl)
  colortbl$y0 = (seq(n)-1)/n
  colortbl$y1 = (seq(n))/n
  
  diff$y0 = colortbl$y0[match(diff$type, names(col_vector) )]
  diff$y1 = colortbl$y1[match(diff$type, names(col_vector) )]
  print(diff)
  # order by color type 
  return(diff)
}

flip_cords <-function(df, paf){
  df = copy(df)
  if(!is.null(paf)){
    major_strand = paf %>% group_by(chr, strand) %>% summarise( total = sum(end-start)) %>% 
      group_by(chr) %>% summarise(flip = strand[which.max(total)])
    df = merge(df, major_strand[, c("chr", "flip")], by="chr", all.x=T)
    if( !("len" %in% colnames(df))){
      df = merge(df, unique(paf[,c("chr", "len")]), by="chr", all.x=T)
    }
  }
  df$tmp_start = df$start
  cond = df$flip == "-"
  df[cond]$start = df$len[cond] - df$end[cond]
  df[cond]$end = df$len[cond] - df$tmp_start[cond]
  
  minus = df$strand == "-"
  plus = df$strand == "+"
  df[cond & minus]$strand = "+"
  df[cond & plus]$strand = "-"
  return(df)
}

hist_data <- function(df, w=1000000){
  ys = c()
  chrs = c()
  starts = c()
  rtn = NULL
  for(chr in unique(df$chr)){
    bed = GRanges( df[df$chr==chr, c("chr","start","end")] )  
    m = max(end(bed))
    
    tmp = data.table(chr=chr, start=seq(1, m, w )); 
    tmp$end = tmp$start + w; 
    tmp$name = paste(paste(tmp$chr, tmp$start, sep=":"), tmp$end, sep="-");
    tmp$range = apply(tmp, 1, function(x) GRanges(x[4]))
    tmp$y = apply(tmp, 1, function(x) sum(width(reduce(pintersect(bed, x[[5]])))) )
    rtn = rbind(rtn, tmp)
    #for(start in seq(1, m, w )){
    #  end = start + w
    #  range = GRanges(c(glue("{chr}:{start}-{end}")))
    #  y = sum(width(reduce(intersect(bed, range))))
    #  ys = c(ys, y); chrs =c(chrs, chr); starts=c(starts, start)
    #}
  }
  #rtn = data.table(chr=chrs, start=starts,y=ys)
  #rtn$end = rtn$start + w
  return(rtn)
}
#x=hist_data(asm_dups[asm_dups$chr %in% c("tig00008002", "tig00018105")]);x

plot_col_data <- function(kp, col){
  if(nrow(col) < 1){ return() }
  m = max(col$coverage)
  kpAxis(kp, r0=0, r1=0.25, tick.pos=round(c(0, 0.5, 1)*m), ymax = m ) 
  kpBars(kp, chr=col$chr, x0=col$start, x1=col$end, y0=0, y1=col$coverage, ymax=m, r0=0, r1=0.25, col="#7f0000", border=NA)
}

plot_dup_hist <- function(kp, dups){
  #kpPlotDensity(kp, data=toGRanges(dups), r0=0.5, r1=1, window.size = 100000)
  w = 1000000
  x = hist_data(dups, w=w);
  kpAxis(kp, r0=0.5, r1=1, tick.pos=round(c(0, 0.5, 1)*w), ymax = w ) 
  kpBars(kp, chr=x$chr, x0=x$start, x1=x$end, y0=0, y1=x$y, ymax=w, r0=0.5, r1=1)

  kpBars(kp, chr=dups$chr, x0=dups$start, x1=dups$end, y0=0, y1=1, r0=0.25, r1=0.5, col=dups$hex, border=NA)
}

plot_diff <- function(kp, diff){
  kpBars(kp, chr=diff$chr, x0=diff$start, x1=diff$end, y0=diff$y0, y1=diff$y1, r0=0, r1=0.25, col=diff$color, border=NA)
}

plot_hg38 <- function(chr, dups=NULL){
  pp = getDefaultPlotParams(plot.type = 5); pp$bottommargin = 0
  kp <- plotKaryotype(chromosomes = chr, genome = "hg38", plot.type = 5, plot.params = pp)
  
  # dups
  if(!is.null(dups)){
    cond = dups$chr == chr
    dups = dups[cond]
    plot_dup_hist(kp, dups)
    kpText(kp, chr=chr, x=-1000000,y=0.5,labels="DupMasker", 
          col = "black", r0=0, r1=.5, cex=1,adj=1)
  }
  
  # gaps
  kpPlotRegions(kp, data=reduce(toGRanges(gaps)), data.panel = "ideogram", col = "#FF6600", r0=-0.05, r1=-.75)
  kpText(kp, chr=chr, x=-1000000,y=0.5,labels="Gaps", 
         data.panel = "ideogram", col = "#FF6600", r0=-0.05, r1=-.75, cex=1,adj=1)
}

plot_asm <- function(chr, paf, dups=NULL, rm=NULL, diff=NULL, col=NULL, plot.type = 3) {


  if(is.null(chr)){
    keep = chrs
  }else{
    keep = unique(paf$chr[paf$prichr == chr ] )
  }
  if(!is.null(rm) && !is.null(paf) ){
    rm = flip_cords( rm[rm$chr %in%  keep ], paf)
  }
  if(!is.null(dups) && !is.null(paf) ){
    dups = flip_cords( dups[dups$chr %in%  keep ], paf)
  }
  if(!is.null(diff) && !is.null(paf) ){
    diff = flip_cords( diff[diff$chr %in%  keep ], paf)
  }
  if(!is.null(col) && !is.null(paf) ){
    col = flip_cords( col[col$chr %in%  keep ], paf)
  }
  if(!is.null(paf) ){
    paf = flip_cords( paf[ paf$chr %in% keep ], paf)
    
    paf$name = paste(paf$t_chr,":",round(paf$t_start/10^6), "-", round(paf$t_end/10^6), sep="")
    
    paf$color = "#9bc2cf"
    paf$color[paf$strand == "-"] = "#e6bbad"
    liftover=toGRanges( paf[,c("chr", "start","end","name","color")] )
    

    # order contigs in karyoplot by there alignemnt to hg38
    chrs = unique(paf[, c("chr", "len")]); chrs$start=0; colnames(chrs) <- c("chr","end","start")
    chrs = chrs[,c("chr", "start","end")]
    chrs = chrs[ chrs$chr %in% keep ]
    
    ordered_chrs = data.table( paf %>% group_by(chr) %>% summarise(ord = min((t_start+t_end)/2)) %>% arrange(ord) )
    ordered_chrs$chr = as.character(ordered_chrs$chr)
    chrs = merge(chrs, ordered_chrs, by="chr", all.x=T)
    chrs = chrs[ order(ord) , ]
    genome <- toGRanges(chrs)
    
    # find sections that do not align to hg38
    unmapped = setdiff(genome, liftover)
    unmapped$name = "na"
    unmapped$gieStain = "#660099"
  }else{
    print("Must have PAFFFFF")
  }
  
  #
  # make karyoplot
  #
  pp = getDefaultPlotParams(plot.type = plot.type); pp$bottommargin = 0; 
  kp <- plotKaryotype(genome = genome, plot.type = plot.type, plot.params = pp, labels.plotter = NULL)
  
  #
  # lift over ideogram 
  #
  if(!is.null(paf) ){
    kpPlotRegions(kp,  data=liftover, data.panel = "ideogram", border = "black", col=liftover$color, avoid.overlapping=FALSE)
    
    # add liftover labels
    kpPlotMarkers(kp, data=liftover, labels = liftover$name, y=.1,  data.panel = 2, cex=0.5,
                  marker.parts= c(0.1,0.4,0.1)) #text.orientation	="horizontal"
    # shown novel bits
    kpPlotRegions(kp, data=reduce(unmapped), data.panel = "ideogram", border = "black", col="#009900")
  }

  #
  # add dup makser to kp
  #
  if(!is.null(dups)){
    plot_dup_hist(kp, dups)
  }
  
  #
  # add nucmer diffs
  #
  if(!is.null(diff)){
    plot_diff(kp, diff)
  }
  
  #
  # add RepeatMasker
  #
  if(!is.null(rm)){
    kpRect(kp, chr = rm$chr, x0=rm$start, x1=rm$end, y0=0, y1=1,  
           data.panel = "ideogram", border=NA, col=rm$color,  r0=-.05, r1=-.75)
  }
  
  #
  # add collapses
  #
  if(!is.null(col) ){
    plot_col_data(kp,col)
  }
  
  # mark new contigs/gaps with orang 
  width = sum(chrs$len)/500
  width = sum(chrs$len)/500
  kpRect(kp, chr = keep, x0=rep(-width, length(keep)), x1=rep(0, length(keep)), y0=0, y1=1, 
         data.panel = "ideogram", border="#FF6600", col="#FF6600",  r0=1.05, r1=1.75)
  
  # add labels
  kpAddChromosomeNames(kp, cex=0.75)
  kpAddBaseNumbers(kp)
} 




# load in stuff
HG38DM = "../Masker/hg38_50k_dupmasker_colors.bed"

RM = "../Masker/chm13.20200611_repeatmasker.out.bed"
DM = "../Masker/chm13.20200611_dupmasker_colors.bed" 
PAF = "../Masker/chm13.20200611_to_hg38.paf"
COL = "../sda_out/coverage/chm13.20200611.collapses.with.cm.bed"


if(F){
  paf = read_paf(PAF)
  asm_dups = read_dupmasker(DM, paf=paf)
  asm_col = fread(COL, col.names = c("chr","start","end","coverage","median","masked","col_len")); 
  asm_rm = read_rm(RM)
#asm_rm = read_rm(RM)
  #hg38_dups = read_dupmasker(HG38DM)
#asm_diff = read_diff(ASMDIFF)
#hg38_diff = read_diff(DIFF)
}


plot_chr <- function(chr, paf, asm_dups=NULL, asm_rm=NULL, asm_col=NULL, hg38_dups=NULL){
  # for some reason args must be global to work with expression
  chr <<- chr; paf <<- paf; 
  asm_dups <<- asm_dups;   asm_rm <<- asm_rm; asm_col <<- asm_col
  hg38_dups <<- hg38_dups
  
  p1=as.ggplot(expression( plot_hg38(chr, dups=hg38_dups) )  ) 
  p2=as.ggplot(expression( plot_asm(chr, paf, dups=asm_dups, rm=asm_rm, col=asm_col) )) 
  g = plot_grid(p1, p2, labels = c('GRCh38', 'Assembly'), label_size = 12, ncol=1, rel_heights = c(1,1.5))
  return(g)
}

plot_chr("chr16", paf, asm_dups = asm_dups, asm_col=NULL, asm_rm=asm_rm, hg38_dups = NULL )
exit()

#return(0)

CHRS = c(paste("chr", seq(1,22),sep="") ,"chrX")


if(T){
  pdf("cmp_asm_hg38_0602.pdf", 16, 9)
  for(chr in CHRS){
    chr <<- chr # for some reason must be global to work
    g = plot_chr(chr, paf, asm_dups = asm_dups, asm_col=asm_col, asm_rm=asm_rm, hg38_dups = hg38_dups )
    print(g)
    print(chr)
  }

  #unused = asm[ !(asm$chr %in% unique(paf$chr)) & asm$len >= 100000]
  #plot_asm(chr=NULL, unused, dups=asm_dups, rm=asm_rm) 
  #plot.new()
  #unused = asm[ !(asm$chr %in% unique(paf$chr)) & asm$len < 100000]
  #plot_asm(chr=NULL, unused, dups=asm_dups, rm=asm_rm) 
  
  dev.off()
  
}else if (F) {
  chr <<- "chr1"
  p1=as.ggplot(expression( plot_hg38(chr, hg38_dups, gaps, hg38_diff) )  ) 
  p2=as.ggplot(expression( plot_asm(chr, asm, paf=paf, dups=asm_dups, rm=asm_rm,  diff=asm_diff) )) 
  g = plot_grid(p1, p2, labels = c('GRCh38', 'Assembly'), label_size = 12, ncol=1, rel_heights = c(1,1.5))
  print(g)
}else{
  chr <<- "chr1"
  paf = read_paf("HiCanu_T2T_strand_seq_to_hg38.paf")
  hg38_diff = read_diff(DIFF)
  asm = read_fai("../strand_seq/saarclust/HiCanu_T2T_simple_merge/clustered_assembly/clustered.fasta.fai")
  gaps = fread(GAPS, col.names = c("chr", "start", "end"))
  
  
  p1=as.ggplot(expression( plot_hg38(chr, hg38_dups, gaps, hg38_diff) )  ) 
  p2=as.ggplot(expression( plot_asm(chr, asm, paf=paf))) 
  g = plot_grid(p1, p2, labels = c('GRCh38', 'Assembly'), label_size = 12, ncol=1, rel_heights = c(1,1.5))
  print(g)
  unused = asm[ !(asm$chr %in% unique(paf$chr)) ]
  print(unused)
}















