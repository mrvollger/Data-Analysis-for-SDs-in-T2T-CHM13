#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")

new_genes = NEW_GENES
new_genes=new_genes[new_genes$cdslen>200,]; dim(new_genes)

textsize=0.8
new_genes = new_genes %>% group_by(`origin gene`) %>% summarise(gene_name=paste(unique(`origin gene`), length(start), sep=" X "), `#chr`=`#chr`[1], start=start[1])
chrs= chrs[chrs %in% unique(new_genes$`#chr`) ]
new_genes$gene_name = gsub(" X 1", "", new_genes$gene_name)
pgenes=NULL
count = 1
for(i in seq(1,length(chrs),3)){
  i2=min(length(chrs), i+2)
  pp <- getDefaultPlotParams(plot.type = 4)
  #pp$topmargin = 0
  #pp$bottommargin = 0
  #pp$data1outmargin = 0
  pgenes[[count]]=as.ggplot(expression(
    kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, chromosomes = rev(rev(chrs)[i:i2]), plot.type = 4),
    kpPlotMarkers(kp, chr = new_genes$`#chr`, x=new_genes$start, labels = new_genes$gene_name, y=0.2, cex=textsize,marker.parts = c(0.1, 0.8, 0.1), 
                  max.iter=1500, ignore.chromosome.ends=T)
  )) #+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  count = count + 1
}
bot = plot_grid(pgenes[[1]], pgenes[[2]], pgenes[[3]], pgenes[[4]],pgenes[[5]],pgenes[[6]], labels=c("c"), nrow=ceiling(length(chrs)/3) ); bot
ggsave("genes.pdf", plot=bot, height = 24/1.5, width = 10/1.5)

ggtitle("New gene copies with ORFs")+theme(plot.title = element_text(hjust = 0.5, face = "bold")); 

pgenes2=as.ggplot(expression(
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, chromosomes = chrs[6:12], plot.type = 4),
  kpPlotMarkers(kp, chr = new_genes$`#chr`, x=new_genes$start, labels = new_genes$gene_name, y=0.2, cex=textsize, marker.parts = c(0.1, 0.8, 0.1), 
                max.iter=1500, label.dist=0.0001,ignore.chromosome.ends	=T)
)
); pgenes2
pgenes3=as.ggplot(expression(
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, chromosomes = chrs[13:length(chrs)], plot.type = 4),
  kpPlotMarkers(kp, chr = new_genes$`#chr`, x=new_genes$start, labels = new_genes$gene_name, y=0.2, cex=textsize, marker.parts = c(0.1, 0.8, 0.1), 
                max.iter=1500, label.dist=0.0001,ignore.chromosome.ends=T	)
)
); pgenes3
