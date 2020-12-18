#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")

inv <- readbed("CHM13_INVs.bed", "INV")
invcol = "darkblue"
Inv=as.ggplot(expression(
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, plot.type = 6),
  kpPlotRegions(kp,NEW, col=transparent(NEWCOLOR,0.5), border = NEWCOLOR),
  kpPlotRegions(kp,reduce(toGRanges(inv)), col=transparent(invcol,0.5), border = invcol),
  legend(x = "bottom", fill = transparent(c(NEWCOLOR, invcol),0.5), legend = c("Unique to CHM13 T2T", "Inversions > 100 kbp"), bty = "n", cex=1.5, ncol=2, inset=0)
)); Inv

ggsave("new_regions_ideo.pdf", height=9, width=16, plot = Inv)