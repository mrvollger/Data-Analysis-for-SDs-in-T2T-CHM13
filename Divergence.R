#!/usr/bin/env Rscript
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")


all = ALL_ALN
#all = all[!grepl('chrY', all$chr)] # filter out chrY since we have no chrY, we now have a y so comment
sds = SDS_ALN
flank = FLANK_ALN
nonsd = NONSD_ALN
chrx = all[all$chr == "chrX"]; chrx$Assembly="chrX"
mhc_rgn = "chr6:28,510,120-33,480,577"
mhc = all[ in_rgn(all, mhc_rgn) ]; mhc$Assembly="MHC"
acro = achro(all); acro$Assembly ="Acrocentric"

df <- bind_rows(sds,flank, nonsd, chrx, mhc,acro); df

pal = c(`SDs`="#af0404", `MHC`="#3282b8", `chrX`="#96bb7c", `SD flanks`="#ede682", `Non SD`="#000000",`Acrocentric`="#ff0000")


df$Assembly = factor(df$Assembly, levels = names(pal) , ordered = TRUE)
df$Divergence = 100-df$perID_by_events
df$aln_frac =  (df$end - df$start) / df$query_length


#
# calculate p values 
#
insync = df[! overlaps(df, NEW)]
p_value = wilcox.test(insync[Assembly == "SDs"]$Divergence,
            insync[Assembly == "Non SD"]$Divergence,
            alternative = "greater")$p.value

if(p_value < 0.0001){
  ptext = c("H[a]:~SDs>Non~SD", 'p<0.0001'); ptext
} else {
  ptext = paste0("H[a]:~SDs>Non~SD:p==",p_value); ptext
} 

fakeadd=0.01
p1 = ggplot() + 
  stat_ecdf(data = insync, 
            aes(Divergence+fakeadd, color = Assembly),
            size=1.5, alpha=0.75) + 
  scale_x_log10(limits=c(fakeadd, NA), breaks = c(fakeadd,0.1,1,10), labels = c("0","0.10","1.00","10.00")) + annotation_logticks(sides="b")+
  scale_fill_manual(values=pal) + scale_color_manual(values=pal) +
  #geom_text(data=NULL,aes(x=1, y=0.25, label = ptext))+
  annotate("text", x=1, y=c(0.25, 0.20),label=ptext, parse=TRUE, hjust=0,size=5)+
  xlab("% divergence of 5kbp windows aligned \n from GRCh38 to T2T CHM13") +
  ylab("Cumulative fraction of 5kbp windows") +
  theme_cowplot()+
  theme(legend.position = "top", legend.title = element_blank()) + guides(fill=guide_legend(ncol=length(pal)/2)) 
p1
ggsave("supp/div_ecdf.pdf", width = 12, height = 8, plot=p1)



div<-toGRanges(df[df$Divergence>=0.5 & df$Assembly=="SDs"])
allsd<-toGRanges(df[df$Assembly=="SDs"])
window=200000
chrs=CHRS
Div=as.ggplot(expression(
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, plot.type = 2, chromosomes=NOYM),
  kpAddBaseNumbers(kp),
  kp2<-kpPlotDensity(kp,div, col=NA,border=NA, window.size = window),
  kpPlotDensity(kp,div, col=transparent(NEWCOLOR,0.25), border=NA,  window.size = window ,ymax=0.1*kp2$latest.plot$computed.values$max.density),
  kp3<-kpPlotDensity(kp,allsd, data.panel = 2, col=NA,border=NA, window.size = window),
  kpPlotDensity(kp, allsd, col=transparent(OLDCOLOR,0.25), border=NA,  window.size = window, ymax=0.1*kp3$latest.plot$computed.values$max.density, data.panel = 2),
  legend(x = "right", fill = transparent(c(NEWCOLOR, OLDCOLOR),0.25), legend = c("SDs > 0.5% diverged", "Density of all SDs"), bty = "n", cex=1, ncol=1, inset=-.0)
))+theme(plot.title = element_text(hjust = 0.5, face = "bold")); Div
ggsave("supp/div_ideo.pdf", width = 12, height = 8, plot=Div)




#
# inversions and new regions compared to GRCh38
#
inv <- readbed("CHM13_INVs.bed", "INV")
invcol = "darkblue"
Inv=as.ggplot(expression(
  trans<-0.2,
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, plot.type = 6, chromosomes = NOYM),
  kpPlotRegions(kp,NEW, col=transparent(NEWCOLOR,trans), border = NEWCOLOR),
  kpPlotRegions(kp,GenomicRanges::reduce(toGRanges(inv)), col=transparent(invcol,trans), border = invcol),
  kpAddCytobands(kp),
  legend(x = "right", 
         fill = transparent(c(NEWCOLOR, invcol, "black"),trans), 
         legend = c("Unique to CHM13 T2T", "Inversions (Strand-seq)", "Centromere"), 
         bty = "n", cex=1, ncol=1, inset=0.05)
)); Inv
ggsave("supp/div_ideo_large_differneces.pdf", width = 12, height = 8, plot=Inv)



#p = plot_grid(p1, Div, bot, Inv,labels = c("a","b",NA,"d"), nrow=2, ncol=2)
p = plot_grid(p1,  
              plot_grid(Div,Inv, labels = c("b","c"), ncol=2),
              labels = c("a", NA), 
              nrow=2);p

ggsave("figures/Divergence.pdf", plot=p, height = 14, width = 14)

      
plot.new = NEW %>%
  mutate(length=end-start)
length_stats(plot.new)
new_lengths = ggplot(data = plot.new,
                     aes(x=length, fill="a")) +
  geom_histogram(bins=25) + 
  scale_x_continuous(label=comma, trans = "log10")+
  annotation_logticks(side="b")+
  ylab("Count") +
  xlab("Genomic length of new or structurally different regions")+
  theme_cowplot() + 
  theme(legend.position = "none")+
  scale_fill_manual(values = c(NEWCOLOR))
xx  = ggdraw(new_lengths)+
  draw_grob(tableGrob(length_stats(plot.new),
                      rows=NULL, theme = ttheme_minimal()),
            x=0, y=.25, width=1, height=1, hjust = 0, vjust = 0)

ggsave("supp/new_region_lengths.pdf", plot=xx, height = 8, width = 12)  







