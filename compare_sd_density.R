setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")

hg_fai = FAI_38 
diff = merge(FAI,hg_fai, by=c("chr"))
diff$scale = diff$chrlen/diff$length
diff[chr=="chr9"]$scale=1
diff[chr=="chr16"]$scale=1


sd_38 = merge(SEDEF_38, diff[, chr, scale]) 
sd_38$start = sd_38$start*sd_38$scale
sd_38$end = sd_38$end*sd_38$scale

con = sd_38$chr=="chr9" & sd_38$start > 50e6 
sd_38[con]$start = sd_38[con]$start  + (diff$chrlen-diff$length)[diff$chr=="chr9"]
sd_38[con]$end = sd_38[con]$end  + (diff$chrlen-diff$length)[diff$chr=="chr9"]


con = sd_38$chr=="chr16" & sd_38$start > 43e6 
sd_38[con]$start = sd_38[con]$start  + (diff$chrlen-diff$length)[diff$chr=="chr16"]
sd_38[con]$end = sd_38[con]$end  + (diff$chrlen-diff$length)[diff$chr=="chr16"]



100*(diff$chrlen - diff$length)/diff$chrlen

chrcounts <<- data.table(
  bind_rows(SEDEF, SEDEF_38) %>%
  group_by(chr, Assembly) %>%
  summarise(count = n(),
            aln_Mbp = sum(end-start) 
            )
)
  # pivot_wider(names_from = "Assembly")
  
z = ggplot(data=chrcounts, aes(fill=Assembly)) + geom_bar(aes(x=Assembly, weight=count )) +
  facet_row(vars(chr)) +
  scale_fill_manual(values = c(NEWCOLOR, OLDCOLOR))+
  theme_cowplot() + 
  ylab(" # of SD")



#
#
# tile to find increased regions
#
#
seqlengths <- FAI$chrlen
names(seqlengths) = FAI$chr
tiles = tileGenome(seqlengths, tilewidth = 5e6, cut.last.tile.in.chrom=T)

count_13 = c()
count_38 = c()
bp_13 = c()
bp_38 = c()
for(i in 1:length(tiles)){
  y = GenomicRanges::findOverlaps(toGRanges(SEDEF), tiles[i])
  z = GenomicRanges::findOverlaps(toGRanges(sd_38), tiles[i])
  count_13 = c(count_13, length(y))
  count_38 = c(count_38, length(z))
  
  bp_13 = c(bp_13, sum(width(GenomicRanges::reduce(toGRanges(SEDEF[queryHits(y)])))))
  bp_38 = c(bp_38, sum(width(GenomicRanges::reduce(toGRanges(sd_38[queryHits(z)])))))
  
  #print(count_38)
  #print(t2t)
}
sd_count = data.table(as.data.frame(tiles))
sd_count$c13 = count_13
sd_count$c38 = count_38
sd_count$bp13 = bp_13
sd_count$bp38 = bp_38
sd_count$diff = sd_count$c13 - sd_count$c38
largest_diff = sd_count %>% arrange(-diff) %>% head(15)

sd_c_hist= ggplot(data=sd_count) +
  geom_histogram(aes(x=log2(c13/c38), fill= diff>=0),
                 breaks=seq(-8,8,0.25), closed="left")+
  theme_cowplot()+ 
  xlab('log2( # CHM13 SDs / # GRCh38 SDs ) per 5Mbp window')+
  theme(legend.position = 'none')+
  scale_fill_manual(values = c(OLDCOLOR, NEWCOLOR))
ggsave(glue('{SUPP}/sd_count_comparison.pdf'), plot = sd_c_hist, height = 8, width = 8)



sd_bp_hist= ggplot(data=sd_count) +
  geom_histogram(aes(x=log2(bp13/bp38), fill= bp13>=bp38),
                 breaks=seq(-8,8,0.25), closed="left")+
  theme_cowplot()+ 
  xlab('log2( CHM13 SD bp / GRCh38 SD bp ) per 5Mbp window')+
  theme(legend.position = 'none')+
  scale_fill_manual(values = c(OLDCOLOR, NEWCOLOR))



#
#
#
#



#/1e6
b=.5
w <<- 1e6*b
ym <<- 1e3*b
pp <<- getDefaultPlotParams(plot.type=3)
pp$ideogramlateralmargin = 0.01
pp$rightmargin <- 0.01
pp$leftmargin <- 0.01
pp$bottommargin <- 0
pp$topmargin <- 0
pp$ideogramheight <-50
pp5 <<- getDefaultPlotParams(plot.type=5)
pp5$rightmargin <- 0.01
pp5$leftmargin <- 0.01
pp5$bottommargin <- 0
pp5$topmargin <- 0
pp5$ideogramheight <-50

sd_density = function(chrs){
  chrs <<- chrs
  p1 =  as.ggplot(expression(
    kp<-plotKaryotype(genome=GENOME, cytobands = CYTOFULL, plot.type = 3, chromosomes = chrs,
                      plot.params = pp),
    kpAddBaseNumbers(kp),
    kpPlotDensity(kp, data=toGRanges(SEDEF), #GenomicRanges::reduce(toGRanges(SEDEF)), 
                   col=transparent(NEWCOLOR, 0.15), 
                  window.size=w, ymax = ym,
                  clipping = F),
    kpPlotDensity(kp, data=toGRanges(sd_38), #GenomicRanges::reduce(toGRanges(sd_38)), 
                   col=transparent(OLDCOLOR, 0.15), 
                  window.size=w, ymax = ym,
                  clipping = F,
                  data.panel = 2, r0=0.05, r1=1.05),
    kpPlotRegions(kp, data=toGRanges(largest_diff), col=transparent('orange', .2), data.panel = "ideogram")
    #kpLines(kp, chr = sd_count$seqnames, x = (sd_count$end+sd_count$start)/2, y=sd_count$diff, ymax = max(sd_count$diff),r0=0, r1=1)
  )) + coord_cartesian(clip = "off")# +theme_nothing() + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
  
  p1
}
dplot = cowplot::plot_grid(NULL,
                  sd_density(NOM[1:4]),
                  sd_density(NOM[5:9]),
                  sd_density(NOM[10:16]),
                  sd_density(NOM[17:24]),
                  ncol=1,
                  rel_heights =c(.25,2,2,2,2)) 


s=0.3
ggsave(glue('{SUPP}/sd_density_comparison.pdf'), plot = cowplot::plot_grid(dplot,
                                                                   cowplot::plot_grid(sd_c_hist, sd_bp_hist, nrow=1, labels = c('b',"c")),
                                                          ncol=1, rel_heights = c(3,1), labels = c('a', NULL)
                                                          ),
       height = s*40, width = s*40)

#+facet_col(vars(seqnames))
