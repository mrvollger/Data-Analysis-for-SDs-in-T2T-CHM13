setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")

seqlengths <- FAI$chrlen
names(seqlengths) = FAI$chr

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
ymax=1#ceiling(log10(max(large_del_ptr$SVLEN)))
labs = c(1)#comma(round(10^seq(bot,ymax)/1e3))

#
# aln data 
#
ptr_cov = fread("data/misc_files/Clint_PTR_mins.50000_maxc.0.1_window.50000_slide.10000.bed.gz", 
                col.names = c("chr", "start", "end"))
ptr_cov$missing = ptr_cov$end - ptr_cov$start
y = has_genes(ptr_cov)
ptr_cov$has_genes = "orange"
ptr_cov$has_genes[-y] = "blue"

ptr_cov = ptr_cov %>% arrange(has_genes, chr, start)
#kp<-plotKaryotype(genome=GENOME, cytobands = CYTOFULL, plot.type = 2)

genes_of_interest = c("ARHGAP11",
                      "NOTCH2",
                      "LPA$",
                      "CYP2D[6,7]",
                      "NOTCH2NLA",
                      "SMN[1,2]",
                      "SRGAP2[B,C]",
                      "TBC1D3[^[:digit:],P]",
                      "BOLA2[^P]*",
                      "TCAF2")
use = c()
for(g in genes_of_interest){
  use=c(use,grep(g, ALL_GENES$gene))
}
plot_genes = data.table(ALL_GENES[use] %>% group_by(gene) %>% summarise(chr=unique(chr), start=min(start), end=max(end)))
plot_genes$human = "black" 
plot_genes$human[overlaps(plot_genes[,2:4], ptr_cov)] = "darkgreen"
plot_genes

sd_density = function(chrs){
  chrs <<- chrs
  p1 =  as.ggplot(expression(
    kp<-plotKaryotype(genome=GENOME, cytobands = CYTOFULL, plot.type = 3, chromosomes = chrs,
                      plot.params = pp),
    kpPlotMarkers(kp, data=toGRanges(plot_genes[,2:4]), labels = plot_genes$gene,
                  label.color = transparent(plot_genes$human, 0.1),
                  ignore.chromosome.ends=T, text.orientation = "vertical", r0=0.05, r1 = 1.5, cex=.85),
    kpPlotRegions(kp, data=toGRanges(ptr_cov), 
                  #y0 = 0, y1 = ptr_cov$sd_diff,
                  #ymin=0, ymax=1,
                  col = transparent(ptr_cov$has_genes, 0.2), 
                  data.panel = 1, r0=0, r1=0.8),
    kpPlotDensity(kp, data=toGRanges(SEDEF[has_genes(SEDEF),]), col=NEWCOLOR,
                  r0=0.1, r1=1.1,
                  data.panel = 2),
    kpAddBaseNumbers(kp),
    kpAddChromosomeNames(kp)
  )) + coord_cartesian(clip = "off")# +theme_nothing() + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
  
  p1
}
dplot = cowplot::plot_grid(NULL,
                           sd_density(NOM[1:4]),
                           sd_density(NOM[5:9]),
                           sd_density(NOM[10:16]),
                           sd_density(NOM[17:24]),
                           ncol=1,
                           rel_heights =c(1,2,2,2,2)) 
s=0.3
ggsave(glue('{SUPP}/nhp_diffs_comparison.pdf'), plot = cowplot::plot_grid(dplot,
                                                                   ncol=1, labels = c('a')
                                                                  ),height = s*40, width = s*40)





#
#
# make a table of maybe human genes
#
#
PTR_MIS = unique(add_genes(ptr_cov, all = TRUE)[, c("chr", "start","end","missing", "gene", "ORF")]) %>% group_by(chr,start,end,missing,gene) %>%
  summarise(`Has ORF`=sum(ORF) > 0);PTR_MIS

names(PTR_MIS)[4] = "Bases missing from Clint_PTR"

PTR_MIS$`has name` = !grepl("^A[L,C][[:digit:]].*", PTR_MIS$gene)
PTR_MIS = PTR_MIS %>% arrange(chr, start, end,-`has name`,gene)

PTR_MIS
