library(ggplot2)
#install.packages('magick')
library(dplyr)
require(gridExtra)
require(scales)
library(RColorBrewer)
library(data.table)
library(ggrepel)
library(glue)
library(cowplot)
library(karyoploteR)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
CHRS= c(paste0("chr",seq(1,22)),"chrX")
#chrom	chromStart	chromEnd	name	score	strandotherChrom	otherStart	otherEnd	otherSize	uid	posBasesHit	testResult	verdict	chits	ccov	alignfile	alignL	indelN	indelS	alignB	matchB	mismatchB	transitionsB	transversionsB	fracMatch	fracMatchIndel	jcK	k2K
#chr1    5694    9627    chr1:248370235  0       _       
#chr1    248370235       248374126       3891    1       1000    
#N/A N/A      N/A     N/A     data/align_both/0018/both0093495        3945    
#13      66      3879    3704    175     119
#56  0.954885279711266        0.95169578622816        0.0465286055796383      0.0467340721401047

wgac_header =c('chr',	'start1',	'end1',	'name',	'score'	,'strand'	,'chr2',	'start2',	'end2',
          'otherSize',	'uid'	,'posBasesHit',	'testResult',	'verdict'	,'chits'	,'ccov'	,'alignfile',	'aln_len'	,'indelN',
          'indelS',	'alnB',	'matchB',	'mismatchB',	'transitionsB',	'transversionsB',	'fracMatch',	'fracMatchIndel',	'jcK',	'k2K')

readwgac=function(f,fai,tag){
  fai = fread(fai, col.names = c("chr","length","x","y","z"))
  df = fread(f, col.names = wgac_header)
  df$Assembly = tag
  df = merge(df, fai, all.x = TRUE)
  df$chr = factor(df$chr, levels =  c(CHRS, unique(df$chr[which(!df$chr %in% CHRS)]) ) , ordered = TRUE)
  df=df[order(chr)]
  return(df)
}
WGAC = readwgac("../WGAC_20200629/data/GenomicSuperDup.tab","../Masker/t2t-chm13.20200727_masked.fasta.fai", tag="WGAC"); WGAC

read = function(f,fai,tag){
  fai = fread(fai, col.names = c("chr","length","x","y","z"))
  df = fread(f, check.names=TRUE);  
  df$Assembly = tag
  colnames(df)[1]="chr"
  
  #df = df[grep("chr8", df$chr)]
  #df = df[df$rm_coverage <= 0.95]
  #if(tag=="T2T CHM13"){
  if(tag!="Celera WGSA"){
    long = fai$chr[fai$length >= 10^6]
    df = df[(df$chr %in% long) & (df$chr2 %in% long)]
  }
  df = merge(df, fai, all.x = TRUE)
  df = df[df$uppercaseMatches >=500]
  
  df$rm_coverage = 100-df$uppercaseMatches/df$aln_len*100
  

  df$chr = factor(df$chr, levels =  c(CHRS, unique(df$chr[which(!df$chr %in% CHRS)]) ) , ordered = TRUE)
  df=df[order(chr)]
  
  return(df)
}

celera = read("../Masker/Celera_WGSA_sedef_out/SDs.browser.bed","../Masker/Celera_WGSA_masked.fasta.fai", "Celera WGSA")
hg38 = read("../Masker/hg38.no_alt_sedef_out/SDs.browser.bed","../Masker/hg38.no_alt_masked.fasta.fai", "GRCh38")
hg37 = read("../Masker/hg19.no_alt_sedef_out//SDs.browser.bed","../Masker/hg19.no_alt_masked.fasta.fai", "GRCh37")
#p733 = read("../Masker/HG00733.hifiasm.trio.hap1_sedef_out/SDs.browser.bed","../Masker/HG00733.hifiasm.trio.hap1_masked.fasta.fai", "HG00733 pat")
#m733 = read("../Masker/HG00733.hifiasm.trio.hap2_sedef_out/SDs.browser.bed","../Masker/HG00733.hifiasm.trio.hap2_masked.fasta.fai", "HG00733 mat")
chm13 = read("../Masker/t2t-chm13.20200727_sedef_out/SDs.browser.bed","../Masker/t2t-chm13.20200727_masked.fasta.fai", "T2T CHM13")

cols <- intersect(colnames(WGAC), colnames(chm13))

df = rbind(hg38[,..cols],hg37[,..cols],chm13[,..cols],celera[,..cols])#,WGAC[,..cols])  #, p733, m733)
colors = c(`T2T CHM13`="#af0404", GRCh38="#3282b8", GRCh37="#414141", `Celera WGSA`="#ede682", `HG00733 pat`="#96bb7c",`HG00733 mat`="#ade498", `WGAC`="#000000")

# remove symetric
print(colnames(df))
cond = ( paste0(df$chr,df$start1) <= paste0(df$chr2,df$start2) ) 
df = df[cond]


pairwise_plot = function(df, facet=FALSE, inter=TRUE){
    if(facet){
      df = df[ (df$chr %in% unique(chm13$chr)) ]
    }
    if(!inter){
      df =df[(df$chr == df$chr2)]
    }
  
    df$cut = cut(100*df$fracMatch,breaks=seq(90,100), include.lowest = TRUE)
    
    p1 = ggplot(data = df) + 
      geom_histogram(aes(100*fracMatch, weight=alnB, fill=Assembly), 
                     alpha=0.9, color="black", binwidth = 0.5, position="dodge") +
      theme_classic() +
      scale_y_continuous(labels = comma)+
      scale_x_continuous( breaks = seq(89.75,99.75,0.5) , labels  = seq(90,100,0.5)) +
      ylab("Sum of aligned bases") + xlab("% identity of pairwise alignments") + 
      guides(fill=guide_legend(ncol=length(colors))) +
      scale_fill_manual(values=colors)
    
    
    p2 = ggplot(data = df) + 
      geom_histogram(aes(alnB, weight=alnB, fill=Assembly), 
                     alpha=0.85, color="black", bins = 20, position="dodge") +
      ylab("Sum of aligned bases") + xlab("Length of pairwise alignments") +
      scale_y_continuous(labels = comma) +
      scale_x_log10(labels = comma) + annotation_logticks(sides="b") +theme(legend.position = "bottom") +
      theme_classic() + theme(legend.position = "none") + scale_fill_manual(values=colors)
    
    if(facet){
      p1 = p1 + facet_wrap(chr ~ ., ncol=3)
      #p2 = p2 + facet_grid(chr ~ .) + theme(legend.position = "top") 
      #p=plot_grid(p1, p2, ncol = 2,labels=c("a","b"))  
      p=plot_grid(p1, ncol = 1, labels=c("a"))
      return(p)
    }
    
    #p3 = ggplot(data = df) + geom_hex(aes(x=length, y=alnB/length)) + 
    #  scale_x_log10() +facet_grid( assembly~.)
    
    p=plot_grid(p1,p2,nrow = 2,labels=c("a","b"))  
    return(p)
}
p=pairwise_plot(df, facet=TRUE, inter=FALSE);p
ggsave("pairwise_bases_intra_by_chr.pdf", plot = p, height = 16, width = 24)

p=pairwise_plot(df, facet=TRUE);p
ggsave("pairwise_bases_all_by_chr.pdf", plot = p, height = 16, width = 24)

p=pairwise_plot(df);p
ggsave("pairwise_bases.pdf",plot = p, height = 12, width = 9)

p=pairwise_plot(df[df$chr == df$chr2]);p
ggsave("pairwise_bases_intra.pdf",plot = p, height = 12, width = 9)

p=pairwise_plot(df[df$chr != df$chr2]);p
ggsave("pairwise_bases_inter.pdf",plot = p, height = 12, width = 9)


p3=ggplot(data = df) +
  geom_histogram(aes(x=rm_coverage,weight=aln_len,color=Assembly,fill=Assembly), binwidth = 1) + 
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  theme_classic()+
  facet_grid(Assembly ~ ., scales="free_y")

#plot_grid(p1,p2,p3, nrow = 3,labels=c("a","b","c"))  





if(F){
contig_lengths = unique(df[,c("chr", "length", "Assembly")])
contig_lengths = contig_lengths[order(Assembly,length)]
contig_lengths$length = as.numeric(contig_lengths$length)
contig_lengths = contig_lengths %>% group_by(Assembly) %>% mutate(cs = cumsum(length))
contig_lengths = contig_lengths %>% group_by(Assembly) %>% mutate(NX = 100-100*cs/sum(length))
contig_lengths
p3 = ggplot(data = contig_lengths) + 
  geom_step(aes(x=NX,y=cs,color=Assembly)) +
  theme_classic() +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) + 
  scale_color_manual(values=colors); 
}

g = unique(df[Assembly=="T2T CHM13",chr,length])
g$end =g$length
g$start=0
genome = toGRanges(data.frame(chr=g$chr, start = g$start, end=g$end))

kp = plotKaryotype(genome = genome)
kpPlotDensity(kp,toGRanges(chm13), window.size=500000)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density)


minid = 0.0
minlen = 0
pdf(glue("SD_links_{minid}_{minlen}.pdf"), height = 12 , width = 16)
kp = plotKaryotype(genome = genome, chromosomes = CHRS)
touse=chm13
inter = touse[fracMatch>minid & aln_len>minlen & chr!=chr2]
intra = touse[fracMatch>minid & aln_len>minlen & chr==chr2]
dim(inter)

kpPlotLinks(kp, 
            data=toGRanges(data.frame(chr=inter$chr,start=inter$start1,inter$end1) ), data2=toGRanges(data.frame(chr=inter$chr2,start=inter$start2,inter$end2)), 
            col=transparent("darkred", amount=0.75), border = NA)
kpAddCytobands(kp)
kpPlotLinks(kp, 
            data=toGRanges(data.frame(chr=intra$chr,start=intra$start1,intra$end1) ), data2=toGRanges(data.frame(chr=intra$chr2,start=intra$start2,intra$end2) ), 
            col=transparent("darkblue", amount=0), border = NA)
kpAddBaseNumbers(kp)

dev.off()
dev.off()
dev.off()


