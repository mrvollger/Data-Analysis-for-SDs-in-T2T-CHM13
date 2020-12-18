#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")

new_genes = read_excel(glue("../Assembly_analysis/Liftoff/{V}.liftoff.summary.xlsx"), sheet="NewCopiesLongestCDS")
new_genes=new_genes[new_genes$cdslen>200,]; dim(new_genes)

all = readbed(glue("../Assembly_analysis/Align/{V}.split.bed"), "All 5kbp windows")
#all = all[!grepl('chrY', all$chr)] # filter out chrY since we have no chrY, we now have a y so comment
sds = readbed(glue("../Assembly_analysis/Align/{V}.split.sd.bed"), "SDs")
flank = readbed(glue("../Assembly_analysis/Align/{V}.split.sdflank.bed"), "SD flanks")
nonsd = readbed(glue("../Assembly_analysis/Align/{V}.split.nosd.bed"), "Non SD")
chrx = all[all$chr == "chrX"]; chrx$Assembly="chrX"
mhc_rgn = "chr6:28,510,120-33,480,577"
mhc = all[ in_rgn(all, mhc_rgn) ]; mhc$Assembly="MHC"
acro = achro(all); acro$Assembly ="Acrocentric"

df <- bind_rows(sds,flank, nonsd, chrx, mhc,acro); df

pal = c(`SDs`="#af0404", `MHC`="#3282b8", `chrX`="#96bb7c", `SD flanks`="#ede682", `Non SD`="#000000",`Acrocentric`="#ff0000")


df$Assembly = factor(df$Assembly, levels = names(pal) , ordered = TRUE)
df$Divergence = 100-df$perID_by_events
df$aln_frac =  (df$end - df$start) / df$query_length



insync = df[! overlaps(df, NEW)]
fakeadd=0.01
p1 = ggplot(data = insync, aes(Divergence+fakeadd, color = Assembly)) + stat_ecdf(size=1.5, alpha=0.75) + 
  scale_x_log10(limits=c(fakeadd, NA), breaks = c(fakeadd,0.1,1,10), labels = c("0","0.10","1.00","10.00")) + annotation_logticks(sides="b")+
  scale_fill_manual(values=pal) + scale_color_manual(values=pal) +
  xlab("% divergence of 5kbp windows aligned from GRCh38 to CHM12 T2T") +
  ylab("Cumulative fraction of 5kbp windows") +
  theme_cowplot()+
  theme(legend.position = "bottom", legend.title = element_blank()) + guides(fill=guide_legend(ncol=length(pal)/2)) ; p1



div<-toGRanges(df[df$Divergence>=0.5 & df$Assembly=="SDs"])
allsd<-toGRanges(df[df$Assembly=="SDs"])
window=200000
chrs=CHRS
Div=as.ggplot(expression(
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, plot.type = 2, chromosomes=NOYM),
  kpAddBaseNumbers(kp),
  kp<-kpPlotDensity(kp,div, col=NA,border=NA, window.size = window),
  kpPlotDensity(kp,div, col=transparent(NEWCOLOR,0.25), border=NA,  window.size = window ,ymax=0.1*kp$latest.plot$computed.values$max.density),
  kp<-kpPlotDensity(kp,allsd, data.panel = 2, col=NA,border=NA, window.size = window),
  kpPlotDensity(kp, allsd, col=transparent(OLDCOLOR,0.25), border=NA,  window.size = window, ymax=0.1*kp$latest.plot$computed.values$max.density, data.panel = 2),
  legend(x = "right", fill = transparent(c(NEWCOLOR, OLDCOLOR),0.25), legend = c("SDs > 0.5% diverged", "Density of all SDs"), bty = "n", cex=1, ncol=1, inset=-.25)
))+theme(plot.title = element_text(hjust = 0.5, face = "bold")); Div
        

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


#
# inversions and new regions compared to GRCh38
#
inv <- readbed("CHM13_INVs.bed", "INV")
invcol = "darkblue"
Inv=as.ggplot(expression(
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, plot.type = 6, chromosomes = NOYM),
  kpPlotRegions(kp,NEW, col=transparent(NEWCOLOR,0.5), border = NEWCOLOR),
  kpPlotRegions(kp,reduce(toGRanges(inv)), col=transparent(invcol,0.5), border = invcol),
  legend(x = "right", fill = transparent(c(NEWCOLOR, invcol),0.5), legend = c("Unique to CHM13 T2T", "Inversions > 100 kbp"), bty = "n", cex=1, ncol=1, inset=0.05)
)); Inv



p = plot_grid(p1, Div, bot, Inv,labels = c("a","b",NA,"d"), nrow=2, ncol=2)
p = plot_grid(plot_grid(Inv,p1,Div,ncol=1,labels = c("a","b","c")), bot, labels = c("a","d"), nrow=1, ncol=2)
ggsave("Divergence.pdf", plot=p, height = 14, width = 14)

      



plot_grid(bot,Div, nol=2)




