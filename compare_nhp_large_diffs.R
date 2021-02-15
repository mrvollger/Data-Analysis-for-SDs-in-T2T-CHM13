setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")

seqlengths <- FAI$chrlen
names(seqlengths) = FAI$chr

#
#
# read in human expanshions
#
#
minsize=1e4
large_del_ptr = fread('../PAV/20210213/results/Clint_PTR/bed/sv_del.bed.gz') %>% filter( END-POS >= minsize )
callable_ptr = fread('../PAV/20210213/results/Clint_PTR/callable/callable_regions_h1_500.bed.gz')

top = max(ceiling(log10(large_del_ptr$SVLEN)))
bot = floor(log10(minsize)) - 1 #min(floor(log10(large_del_ptr$SVLEN)))
large_del_ptr$color = c('darkblue', "orange", NEWCOLOR)[ceiling(log10(large_del_ptr$SVLEN))-(bot+1)]

large_del_38 = fread('../PAV/20210213/results/GRCh38_chrOnly//bed/sv_del.bed.gz') %>% filter( END-POS >= minsize )
callable_38 = fread('../PAV/20210213/results/GRCh38_chrOnly//callable/callable_regions_h1_500.bed.gz')
large_del_38$color = c('darkblue', "orange", NEWCOLOR)[ceiling(log10(large_del_38$SVLEN))-(bot+1)]



#/1e6
b=.5
w <<- 1e6*b
ym <<- 1e3*b
pp <<- getDefaultPlotParams(plot.type=3)
pp$ideogramlateralmargin = 0.01
pp$rightmargin <- 0.01
pp$leftmargin <- 0.05
pp$bottommargin <- 0
pp$topmargin <- 0
pp$ideogramheight <-50
pp5 <<- getDefaultPlotParams(plot.type=5)
pp5$rightmargin <- 0.01
pp5$leftmargin <- 0.01
pp5$bottommargin <- 0
pp5$topmargin <- 0
pp5$ideogramheight <-50
ymax=ceiling(log10(max(large_del_ptr$SVLEN)))
labs = comma(round(10^seq(bot,ymax)/1e3))
sd_density = function(chrs){
  chrs <<- chrs
  p1 =  as.ggplot(expression(
    kp<-plotKaryotype(genome=GENOME, cytobands = CYTOFULL, plot.type = 3, chromosomes = chrs,
                      plot.params = pp),

    kpRect(kp, data=toGRanges(large_del_38), 
                  y0=bot, y1=log10(large_del_38$SVLEN),
                  ymin=bot, ymax=ymax,
                  col= large_del_38$color,
                  border = large_del_38$color,
                  #window.size=w,
                  clipping = F,
                  data.panel = 1, r0=0.1, r1=1.1),
    kpAxis(kp, ymin=bot, ymax=ymax, labels = labs, tick.pos = seq(bot,ymax),
           data.panel = 1, r0=0.1, r1=1.1),
    kpPlotRegions(kp,  data=GenomicRanges::gaps(toGRanges(callable_38),start=1L, end=seqlengths(toGRanges(callable_38))),
                  col=transparent(NEWCOLOR, .2),
                  data.panel = 1, r0=0.05, r1=0.1),
    
    #
    # PTR
    #
    kpRect(kp, data=toGRanges(large_del_ptr), 
           y0=bot-1, y1=log10(large_del_ptr$SVLEN),
           ymin=bot-1, ymax=ymax,
           col= large_del_ptr$color,
           border = large_del_ptr$color,           #window.size=w,
           clipping = F,
           data.panel = 2, r0=0.1, r1=1.1),
    kpAxis(kp, ymin=bot, ymax=ymax, labels = labs, tick.pos = seq(bot,ymax),
           data.panel = 2, r0=0.1, r1=1.1),
    kpPlotRegions(kp, data=GenomicRanges::gaps(toGRanges(callable_ptr),start=1L, end=seqlengths(toGRanges(callable_ptr))),
                  col=transparent(NEWCOLOR, .2),
                  data.panel = 2, r0=0.05, r1=0.1),

    #
    #
    #
    kpAddBaseNumbers(kp)
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
ggsave('supp/nhp_diffs_comparison.pdf', plot = cowplot::plot_grid(dplot,
                                                                   ncol=1, labels = c('a')
                                                                  ),
       height = s*40, width = s*40)

#+facet_col(vars(seqnames))
