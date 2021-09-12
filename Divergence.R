#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(circlize)
#source("plotutils.R")
pal = c(`MHC`="#3282b8",
        `chrX`="#96bb7c",
        `Mixed` = "#708090",
        `Non SD`="#000000",
        `SDs`="#af0404",
        `SD flanks`="#ede682",
        `Acrocentric`="#ff0000"
       )

WINDOW_SIZE = 10e3
STEP=1e3

#
# general regions
#
notSD = gaps(NONR_SD, start=1L, end=seqlengths(NONR_SD))
notNEW = gaps(GRanges(NEW), start=1L, end=seqlengths(GRanges(NEW)))
MHC_GR = GRanges(data.table(chr="chr6", start=28510120, end=33480577))
CHRX_GR = GRanges(data.frame(chr="chrX", start=0, end=FAI[chr=="chrX"]$chrlen))
Region_Size = data.table(Region = c("Non SD", "SDs", "chrX", "MHC", "Mixed"),
                   Rsize = c(bedlength(intersect(notSD, notNEW)),
                             bedlength(intersect(NONR_SD, notNEW)),
                             bedlength(intersect(CHRX_GR, notNEW)),
                             bedlength(intersect(MHC_GR, notNEW)),
                             0
                   )
)


#
# read in all PAV data
#
allpav = copy(ALLPAV) %>% filter(chr != "chrY") %>% data.table()
allpav$Region = allpav$Assembly


#
# make the windows over the genome to summarize pav results 
#
summarise_pav_by_window <- function(in_df){
  #gr.windows <- tileGenome(Seqinfo(seqnames = as.character(FAI$chr), seqlengths = FAI$chrlen), 
  #                         tilewidth=WINDOW_SIZE, cut.last.tile.in.chrom=TRUE)
  gr.windows <- unlist(
    slidingWindows(
      GRanges(data.table(seqnames = as.character(FAI$chr), 
                         start = 0, 
                         end = FAI$chrlen)), 
      width = WINDOW_SIZE, step = STEP)
    )
  
  # filter out chrY and chrM
  gr.windows = gr.windows[seqnames(gr.windows) %in% NOYM]
  window.df = as.data.table(gr.windows); window.df$windowID = seq(nrow(window.df))
  
  # find overlaps with PAV resutls, consider only the starting base of event
  pav_start = in_df[,c(1,2,2)]; colnames(pav_start)[3] ="end"
  fo <- findOverlaps(GRanges(pav_start), gr.windows, select="all", type="any"); fo
  
  pav_hits = data.table(pavID = queryHits(fo), windowID = subjectHits(fo))
  pav_hits[, colnames(in_df)] = in_df[pav_hits$pavID]
  pav_w = merge(window.df, pav_hits, by="windowID", all.x=T) %>% arrange(windowID) %>% data.table()
  
  pav_windows = pav_w %>% 
    group_by(seqnames, start.x, end.x, windowID) %>%
    summarise(`# SNVs` = sum(SVTYPE == "SNV", na.rm = T),
              `# INS` = sum(SVTYPE == "INS", na.rm = T),
              `# DEL` = sum(SVTYPE == "DEL", na.rm = T),
              `# INV` = sum(SVTYPE == "INV", na.rm = T)
              ) %>%
    data.table(); colnames(pav_windows)[1:3]=c("chr", "start","end")
  
  pav_windows$Divergence = 100*rowSums(pav_windows[,c("# SNVs", "# INS", "# DEL", "# INV")])/(pav_windows$end - pav_windows$start + 1)
  pav_windows
}

add_regions_to_data <- function(rgns){
  rgns$Region = "Mixed"
  rgns[ overlaps(rgns, as.data.table( notSD ), mincov=.90), ]$Region="Non SD"
  rgns[ overlaps(rgns, as.data.table(NONR_SD), mincov=.90), ]$Region = "SDs"
  
  chrx = rgns[rgns$chr == "chrX"]; chrx$Region="chrX"
  
  mhc_rgn = "chr6:28,510,120-33,480,577"
  mhc = rgns[ in_rgn(rgns, mhc_rgn) ]; mhc$Region="MHC"
  
  rgns = rbind(rgns, chrx, mhc)
  rgns$Region = factor(rgns$Region, levels = c("Mixed", "Non SD", "SDs", "chrX", "MHC"))
  rgns
}

#
# windowed dataframe
#
pav_windows = summarise_pav_by_window(allpav)
div.df = add_regions_to_data(pav_windows)
div.df$non_overlap_end = div.df$start + STEP
insync = div.df[!overlaps(div.df, NEW)]

out.tbl = insync[order(as.character(chr), start),][
  Region %in% c("SDs", "Non SD"), c("chr", "start", "non_overlap_end", "# SNVs")]
fwrite(file = "data/misc_files/SNV.pav.bed.gz", out.tbl,
            col.names = F, row.names = F, sep="\t", quote = F, scipen = 20)
#
# all data 
#
pav = add_regions_to_data( allpav[overlaps(allpav, as.data.table(notNEW))] )
pav_start = pav[,c(1,2,2)]; colnames(pav_start)[3] ="end"
fo = findOverlaps(GRanges(pav_start)+1, notNEW)
pav$SyntenicRegion = "None"
pav[queryHits(fo)]$SyntenicRegion = as.character( notNEW[ subjectHits(fo) ] )


divergence_summary = pav %>% 
  filter(Region != "Mixed") %>%
  left_join(Region_Size) %>%
  group_by(Region) %>% 
  summarise(`# SNV` = sum(SVTYPE=="SNV"),
            `# INS` = sum(SVTYPE=="INS"),
            `# DEL` = sum(SVTYPE=="DEL"),
            Mbp = unique(Rsize)/1e6,
            `SNVs per kbp` = `# SNV` / (Mbp * 1000), #1000 * median(mismatches / (end - start)), #
            `INSs per kbp` = `# INS` / (Mbp * 1000),
            `DELs per kbp` = `# DEL` / (Mbp * 1000)
  ) %>% data.table(); divergence_summary
table(pav[SVTYPE=="SNV"]$Region)

#
# filter the pav table by snvs
#


div.df[which.max(`# SNVs`)]
insync[which.max(`# SNVs`)]
bedlength(notNEW)/1e6
bedlength(insync[Region=="SDs"])/1e6
bedlength(insync[Region=="Non SD"])/1e6
bedlength(insync[Region=="Mixed"])/1e6


p_greater = wilcox.test(insync[Region == "SDs"]$Divergence,
                      insync[Region == "Non SD"]$Divergence,
                      alternative = "greater")


p_value = wilcox.test(insync[Region == "SDs"]$Divergence,
            insync[Region == "Non SD"]$Divergence)$p.value

p_value = p_greater$p.value
if(is.na(p_value) || p_value < 0.0001){
  ptext = c("H[a]:Non~SD < ~SDs", 'p[value]<0.0001'); ptext
} else {
  ptext = c("H[a]:Non~SD < ~SDs",  paste("p[value] ==",signif(p_value,2)) ); ptext
} 

fakeadd=0.01
p1 = ggplot() + 
  stat_ecdf(data = insync %>% filter(Region != "Mixed"), 
            aes(Divergence+fakeadd, color = Region),
            size=1.5, alpha=0.75) + 
  scale_x_log10(limits=c(fakeadd, 100), 
                breaks = c(fakeadd,0.1,1,10), 
                labels = c("0","0.10","1.00","10.00")
                ) +
  annotation_logticks(sides="b")+
  scale_fill_manual(values=pal) + scale_color_manual(values=pal) +
  annotate("text", x=1, y=c(0.25, 0.2),label=ptext, parse=TRUE, hjust=0,size=5)+
  xlab(glue("% divergence of {WINDOW_SIZE/1e3}kbp windows ({STEP/1e3}kbp slide) aligned \n from GRCh38 to T2T CHM13")) +
  ylab(glue("Cumulative fraction of {WINDOW_SIZE/1e3}kbp windows")) +
  theme_cowplot()+
  theme(legend.position = "top", legend.title = element_blank()) + 
  guides(fill=guide_legend(ncol=length(pal)/2)) 
p1
ggsave(glue("{SUPP}/div_ecdf.pdf"), width = 12, height = 8, plot=p1)



sd_vs_non_p_hist = ggplot(data=insync %>%
         filter(Region %in% c("SDs", "Non SD")) 
       ) + 
  geom_histogram(aes(x=`# SNVs`+0.1, fill=Region, y=..density..), 
                 bins = 60, position = "dodge")+
  #, binwidth = 1) +
  #geom_density_ridges(aes(x=`# SNVs`+0.1, fill=Region, y = Region))+
  #geom_density(aes(x=snps_per_kbp+0.0, fill=Region, alpha=0.5))+
  #geom_text(data = . %>% group_by(Region) %>%
  #            summarise(mean = round(mean(snps_per_kbp),2), median =round(median(snps_per_kbp),3)),
  #          aes(x= mean, y= 2, label = paste(mean, median, sep=", "), color=Region)
  #          ) +
  scale_x_log10(limits=c(0.09, NA), 
                breaks = c(0.1,1,10,100), 
                labels = c("0","1","10", "100")
                ) + annotation_logticks(sides="b") +
  #scale_y_log10(labels = comma) +
  #facet_col(~Region, scales = "free_y") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
    #coord_cartesian(xlim = c(0,10)) +
  xlab(glue("# of SNVs per {WINDOW_SIZE/1e3}kbp window ({STEP/1e3}kbp slide) aligned \n from GRCh38 to T2T CHM13")) +
  theme_cowplot() + 
  theme(legend.position = "top"); sd_vs_non_p_hist
ggsave(glue("{SUPP}/SNVs_per_kbp_SD_vs_Not.pdf"), plot=sd_vs_non_p_hist, height = 8, width = 12)


#
# DISTANCE between SNVs
#
SNVs = pav %>% 
  group_by(chr, SyntenicRegion) %>%
  filter(chr %in% NOYM & SVTYPE == "SNV" & ! Region %in% c("Mixed", "MHC", "chrX") ) %>%
  mutate( dist = start - lag(start)) %>% 
  #filter(dist < 1e6) %>%
  data.table(); SNVs
table(SNVs$Region)

SNV_dist = ggplot(data = SNVs) + 
  geom_vline(data=data.frame(), aes(xintercept=seq(1,10)), size=0.2, linetype = "dashed") +
  #geom_histogram(aes(x=dist, color=Region, fill=Region), #, y=..density..), 
  #               bins = 50)+ #, position = "dodge") +
  #geom_density_ridges(aes(x=dist, color=Region, y=Region,  fill=Region), alpha=0.85, bandwidth = 0.015) +
  geom_density(aes(x=dist, y=..scaled.., color=Region,  fill=Region), alpha=0.4, adjust = .5, size=1) +
  geom_rug(data = . %>% filter(Region=="SDs"), aes(x=dist, color=Region), alpha=0.4)+
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  #facet_zoom(x, xlim = c(3,5), shrink = F)+
  scale_x_log10(breaks = 10^seq(0,5),
                labels = comma) +  annotation_logticks(side="b") +
  #scale_y_log10(labels = comma) +
  #facet_col(~Region, scales="free_y")+
  xlab("Distance between consecutive SNVs")+
  theme_cowplot() + 
  theme(legend.position = "top"); SNV_dist
ggsave(glue("{SUPP}/SNV_distance.pdf"), plot=SNV_dist, height = 8, width = 12)



#
# SNVs ideogram
#
my_ideo_snv <<- function(){
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTOFULL, plot.type = 1, chromosomes=NOYM)#c("chr1"))
  #kpAddLabels(kp, labels = "Number SNV per kbp", label.margin = 0.06, srt=90, pos=3)
  nparts <- c(0,2,5,Inf)
  heatcol = colorRampPalette(c(OLDCOLOR, NEWCOLOR))(4)
  for(i in 1:(length(nparts)-1)) {
    at <- autotrack(i, length(nparts))
    kpAddLabels(kp, labels = paste(nparts[i],nparts[i+1], sep = "-"), r0=at$r0, r1=at$r1, cex=0.5)
    #kpPlotHorizon(kp, data=rand.data, num.parts = nparts[i],
    #   ymin=ymin, ymax=ymax,
    kpPlotRegions(kp, data = toGRanges(insync[Region == "SDs" & `# SNVs` > nparts[i] & `# SNVs` <= nparts[i+1]] ), 
                  r0=at$r0, r1=at$r1, avoid.overlapping = F,
                  col= heatcol[i])
  }
}
snvs_per_kbp = as.ggplot(expression(my_ideo_snv()))
ggsave(glue("{SUPP}/snvs_per_kbp.pdf"), plot=snvs_per_kbp, height = 16, width = 12)




#
#
#
#
div<-toGRanges(insync[insync$Divergence>=0.2 & insync$Region=="SDs"])
allsd<-toGRanges(insync[insync$Region=="SDs"])
window=200000
chrs=CHRS
Div=as.ggplot(expression(
  kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, plot.type = 2, chromosomes=NOYM),
  kpAddBaseNumbers(kp),
  kp2<-kpPlotDensity(kp, div, col=NA,border=NA, window.size = window),
  kpPlotDensity(kp,div, col=transparent(NEWCOLOR,0.25), border=NA,  window.size = window), #ymax=0.1*kp2$latest.plot$computed.values$max.density),
  kp3<-kpPlotDensity(kp,allsd, data.panel = 2, col=NA,border=NA, window.size = window),
  kpPlotDensity(kp, allsd, col=transparent(OLDCOLOR,0.25), border=NA, data.panel = 2, window.size = window),# ymax=0.1*kp3$latest.plot$computed.values$max.density, data.panel = 2),
  legend(x = "right", fill = transparent(c(NEWCOLOR, OLDCOLOR),0.25), legend = c("SDs > 0.2% diverged", "Density of all SDs"), bty = "n", cex=1, ncol=1, inset=-.0)
))+theme(plot.title = element_text(hjust = 0.5, face = "bold")); Div
ggsave(glue("{SUPP}/div_ideo.pdf"), width = 12, height = 8, plot=Div)



#
#
# SNV heatmap as circous
#
#
merged = as.data.table(reduce(GRanges(insync[Region == "SDs"])+1e5)) %>% 
  mutate(sd_id = paste(seqnames,start,end,sep="_"),
         dist_bewteen = start - lag(end),
         sd_start=start,
         sd_end = end) %>% 
  mutate(start)
merged[start<0, c("start","sd_start")] = 0
merged
bedlength(merged)/1e6
o = findOverlaps(GRanges(insync[Region %in% c("SDs", "Non SD")]), GRanges(merged), minoverlap = 1e4)
heatmap.df = cbind(insync[Region %in% c("SDs", "Non SD")][queryHits(o)]
                   ,merged[subjectHits(o), c("sd_id", "sd_start", "sd_end")])
heatmap.df$color = NEWCOLOR
heatmap.df[Region != "SDs"]$color = OLDCOLOR

f1 = function(){
  cyto.df <- data.frame(CYTOFULL)[, c("seqnames","start","end","seqnames","gieStain")]
  cyto.df = cyto.df[cyto.df$seqnames %in% touse,]
  cyto.df$seqnames = factor(cyto.df$seqnames, levels =  c(CHRS[!CHRS %in% ACHRO ], ACHRO))
  cyto.df  = cyto.df[order(cyto.df$gieStain, cyto.df$seqnames),]
  cyto.df$seqnames = as.character(cyto.df$seqnames)
  #gap.degree=360/(4*length(unique(seqnames(CYTOFULL)))); gap.degree
  #circos.par(cell.padding = c(0, 0, 0, 0), gap.degree=gap.degree)
  circos.initializeWithIdeogram(cyto.df, ideogram.height = .1, plotType = c("ideogram", "labels"))
  
  col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))
  bed=insync[chr %in% touse, c("chr","start","end","# SNVs")]
  bed
  circos.genomicHeatmap(bed, col = col_fun, side = "inside", border = "white")
  
};
f2 = function(){
  conversion = conversion[chr %in% touse]
  print(dim(conversion))
  gap.degree=360/(8*length(conversion$sd_id)); gap.degree
  circos.par(cell.padding = c(0, 0, 0, 0), track.height=0.4, gap.degree=gap.degree)
  circos.initialize(factors=conversion$sd_id, xlim = conversion[,c(2,3)])
  zoom_data = heatmap.df[sd_id %in% conversion$sd_id][`# SNVs` > 10]
               
  circos.track(zoom_data$sd_id, x = zoom_data$start, y = log10(zoom_data$`# SNVs`+1), 
               panel.fun = function(x, y, z) {
                 circos.points(x, y, pch = 16, cex = 0.5, col=NEWCOLOR)
               },
               track.margin = c(0, 0)
               )
  
  o = findOverlaps(reduce(GRanges(SEDEF[chr %in% touse])), GRanges(merged))
  SD = cbind(as.data.table(reduce(GRanges(SEDEF[chr %in% touse])))[queryHits(o)], 
             merged[subjectHits(o), c("sd_id", "sd_start", "sd_end")])
  SD[start < sd_start]$start = SD[start < sd_start]$sd_start
  SD[end > sd_end]$end = SD[end > sd_end]$sd_end
  circos.par(cell.padding = c(0, 0, 0, 0), track.height=0.1)
  circos.track(ylim = c(0,1))
  circos.genomicRect(SD[,c("sd_id", "start", "end")], 0, col=transparent(NEWCOLOR,0.5), border = NA, track.height=0.1)
  
  
  #o = findOverlaps(GRanges(ALLPAV[chr %in% touse]), GRanges(merged))
  #snvs = cbind(ALLPAV[chr %in% touse ][queryHits(o)]
  #           ,merged[subjectHits(o), c("sd_id", "sd_start", "sd_end")])
  #circos.par(cell.padding = c(0, 0, 0, 0), track.height=0.1)
  #circos.track(ylim = c(0,1))
  #circos.genomicDensity(as.data.table(zoom_data[,c("sd_id","start","end")]))
  
  dim(SD)
} ; f2()

conversion=unique(heatmap.df[,c("sd_id","sd_start","sd_end","chr")])
conversion$start2=conversion$sd_start
conversion$end2 = conversion$sd_end
conversion

touse <<- c("chr1")
circos.nested(f2, f1, conversion, connection_height = mm_h(20),connection_col = "lightgray")



#### PDF booklet of SNV density
NONR_SD 
sd_regions = as.data.table( reduce(NONR_SD[width(NONR_SD) > 5e4] + 1e5)  )[width > 5e5 & seqnames %in% NOYM]
colnames(sd_regions) = c("sdb_chr", "sdb_st", "sdb_en", "sdv_width", "sdb_strand")
sd_regions[sdb_st < 0,]$sdb_st = 0
fwrite(sd_regions[,1:3], file="data/misc_files/large_sd_regions.bed.gz", quote = F,
       col.names=F, row.names=F, sep="\t", scipen = 1000)
dim(sd_regions); #bedlength(sd_regions)/1e6
pav_marked_by_block = merged_overlaps(pav, sd_regions) %>% 
  separate(SyntenicRegion, 
           into = c("syn_chr","syn_st","syn_en"), 
           sep = ":|-", convert = T) %>% 
  data.table()

pav_marked_by_block$tag = paste(pav_marked_by_block$sdb_chr, 
                                pav_marked_by_block$sdb_st, 
                                pav_marked_by_block$sdb_en,
                                sep = "_")
pav_snv_plot_df = pav_marked_by_block[tag == pav_marked_by_block$tag[1000]]
ggplot(data=pav_snv_plot_df %>% filter(SVTYPE=="SNV")) +
  geom_histogram(aes(x = start, fill=Region), binwidth = 1e4) +
  scale_color_manual(values = c(Syntenic = "darkgreen")) +
  geom_segment(aes(y = -1, yend=-1,
                   x=syn_st, xend=syn_en, color="Syntenic"), alpha=0.8, size=4) +
  scale_x_continuous(label=comma)+
  xlab("genomic position") +
  ylab("# SNVs in 10 kbp widnows") +
  theme_cowplot() +
  ggtitle(pav_snv_plot_df$tag[1])+
  coord_cartesian(ylim=c(-0.25,NA), 
                  xlim=c(min(pav_snv_plot_df$sdb_st),
                         max(pav_snv_plot_df$sdb_en))
                  ) 

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
))
ggsave(glue("{SUPP}/div_ideo_large_differneces.pdf"), width = 12, height = 8, plot=Inv)



#p = plot_grid(p1, Div, bot, Inv,labels = c("a","b",NA,"d"), nrow=2, ncol=2)
p = cowplot::plot_grid(
              cowplot::plot_grid(p1, sd_vs_non_p_hist, SNV_dist, labels = c("a","b","c"), ncol=1),  
              cowplot::plot_grid(snvs_per_kbp, labels = c("d")),
              ncol=2)
ggsave(glue("{SUPP}/Divergence.pdf"), plot=p, height = 24/1.5, width = 20/1.5)


z1 = ggplot(data=pav %>% filter(chr == "chr16" & start > 20e6 & end < 60e6))+
  geom_point(aes(x=start, y = runif(length(start)), color=Region) ) + 
  theme_cowplot() + 
  theme(legend.position = "top"); z1

z2 = ggplot(data=pav %>% filter(SVTYPE == "SNV" & chr == "chr16" & start > 15e6 & end < 20e6))+
  geom_histogram(aes(x=start, color=Region, fill = Region), binwidth = 10000 ) + 
  theme_cowplot() +
  #scale_y_log10()+
  theme(legend.position = "top")
z2

write.table(pav[SVTYPE=="SNV" & Region %in% c("SD", "Non SD"), 1:3],
          file = "~/Desktop/pav.snv.bed",
          sep="\t", quote = F, row.names = F, col.names = F)


#
# length of non syntenic regions
#
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

ggsave(glue("{SUPP}/new_region_lengths.pdf"), plot=xx, height = 8, width = 12)  











##########################################
##########################################
##########################################
##########################################
########### OLD BAD CODE #################
##########################################
##########################################
##########################################
if(F){
    all = ALL_ALN
    all = all[!grepl('chrY', all$chr)] # filter out chrY since we have no chrY, we now have a y so comment
    sds = ALL_ALN[overlaps(ALL_ALN, as.data.table(NONR_SD), mincov=.95)]; sds$Assembly="SDs" #SDS_ALN
    flank = FLANK_ALN
    nonsd = NONSD_ALN
    chrx = all[all$chr == "chrX"]; chrx$Assembly="chrX"
    mhc_rgn = "chr6:28,510,120-33,480,577"
    mhc = all[ in_rgn(all, mhc_rgn) ]; mhc$Assembly="MHC"
    acro = achro(all); acro$Assembly ="Acrocentric"
    
    df <- bind_rows(sds,flank, nonsd, chrx, mhc,acro); df
    
    
    
    df$Assembly = factor(df$Assembly, levels = names(pal) , ordered = TRUE)
    df$Divergence = 100-df$perID_by_events
    df$aln_frac =  (df$end - df$start) / df$query_length
  
  #
  # calculate p values 
  #
  insync = df[! overlaps(df, NEW)] %>%
    filter((end-start) > 2000) %>% data.table()
  
  
  #
  #
  # SNP calculations 
  #
  #
  bed <- GRanges(readbed("data/misc_files/tmp.bam.bed.gz", "z"))
  cov = coverage(bed )
  
  
  #ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
  
  
  
  tmpSNVs = readbed("data/misc_files/chm13.draft_v1.0_plus38Y.snv.bed.gz", "Non SD")
  tmpSNVs = readbed("data/misc_files/snv_snv.bed.gz", "Non SD")
  tmpSNVs = readbed("data/misc_files/snv_snv.bed.gz", "Non SD")
  tmpSNVs = readbed("data/misc_files/snv_snv.bed.gz", "Non SD")
  
  
  tmpSNVs = tmpSNVs[! overlaps(tmpSNVs, NEW)] 
  tmpSNVs[overlaps(tmpSNVs, as.data.table(NONR_SD))]$Assembly = "SDs"
  
  SNVs = tmpSNVs %>% 
    mutate(Region = Assembly) %>%
    filter(chr %in% NOYM) %>%
    group_by(chr) %>%
    mutate( dist = end - lag(start)) %>% 
    #separate(V9, c("DP", "info"), sep=";") %>%
    #filter( DP == "DP=1") %>%
    filter(dist < 1e6) %>%
    data.table(); SNVs
  
  #
  #
  # make tables and stuff for evans summary statistics 
  #
  insync = df[! overlaps(df, NEW)] %>% mutate(snps_per_kbp = 1000 * mismatches / (end-start),
                                              ins_per_kbp = 1000 * insertion_events / (end-start),
                                              del_per_kbp = 1000 * deletion_events / (end-start)) %>%
    arrange(chr,start,end) 
  insync$Region = insync$Assembly
  #
  #
  #
  window_averages = function(zz, window_size, colname){
    gr.windows <- tileGenome(Seqinfo(seqnames = as.character(FAI$chr), seqlengths = FAI$chrlen), 
                             tilewidth=window_size, cut.last.tile.in.chrom=TRUE)
    gr.data <- GRanges(zz)
    gr.data.cov <- GenomicRanges::coverage(gr.data, weight=colname) 
    seqlevels(gr.windows, pruning.mode="coarse") <- names(gr.data.cov)
    
    # windows with actual coverage 
    gr.windows.with.cov = intersect(gr.windows, GenomicRanges::reduce(gr.data+1000) )
    
    binnedSum(gr.windows, gr.data.cov, colname)
  }
  #sd_windows = window_averages(insync[Assembly=="SDs"], 1e5, "snps_per_kbp")
  #non_windows = window_averages(insync[Assembly=="Non SD"], 1e5, "snps_per_kbp")
  
  
  divergence_summary = insync %>% 
    filter(Region != "Mixed") %>%
    left_join(Region_Size) %>%
    group_by(Region) %>% 
    summarise(`# SNV` = sum(`# SNVs`),
              `# INS` = sum(`# INS`),
              `# DEL` = sum(`# DEL`),
              Mbp = unique(Rsize)/1e6,
              `SNVs per kbp` = `# SNV` / (Mbp * 1000), #1000 * median(mismatches / (end - start)), #
              `INSs per kbp` = `# INS` / (Mbp * 1000),
              `DELs per kbp` = `# DEL` / (Mbp * 1000)
    ) %>% data.table(); divergence_summary
  
}
