#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")
#sedef=readbed("../Assembly_analysis/SEDEF/{V}.SDs.bed", "T2T CHM13")
#duplicons = readbed("../Assembly_analysis/Masked/{V}_dupmasker_colors.bed", "T2T CHM13")
#enriched = readbed("../Assembly_analysis/SEDEF/{V}.sedef.enriched.bed", "highsd")[,1:4]
#rm = readbed("../Assembly_analysis/Masked/{V}_repeatmasker.out.bed", "T2T CHM13", rm=T)
rm=RM
duplicons=DM_BED
sedef=SEDEF
probes = readbed("data/misc_files/probe.mappings.bed.gz", "probes")
probes = probes[!probe_id %in% c("174552_ABC10_2_1_000044707300_L19::chr22:11387009-11432887", "174222_ABC10_2_1_000044616600_A12::chr14:7053457-7120716")]
pnames  = data.table(probe_id=unique(probes$probe_id))
pnames$color = brewer.pal(length(pnames$probe_id), "Dark2")
probes = data.table(merge(probes, pnames, by ="probe_id", all.x=T))
probes = data.table(probes %>% relocate(probe_id, .after=end))
probes$label = as.numeric(as.factor(probes$probe_id))
colors = pnames$color
names(colors)=pnames$probe_id
sat = rm[rm$type == "Satellite"]


rm=rm[end-start>10000]
dim(rm)


enriched_sds=data.frame(reduce(toGRanges(probes)+5e6))
colnames(enriched_sds)[1]="chr"
enriched_sds[enriched_sds$start < 1,]$start = 1
names=unique(enriched_sds$chr)
duplicons=duplicons[chr %in% names]
duplicons$hex = sapply(strsplit(duplicons$V9, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))


zoom_rgns  = data.table(enriched_sds %>% 
  group_by(chr) %>% 
  #summarise(start=min(start), end=max(end)) %>% 
  mutate(Mbp=(end-start)/10^6, name=paste0(chr,": ",round(start/10^6),"-",round(end/10^6), " Mbp")) )
zoom_rgns$chr = factor(zoom_rgns$chr, levels =  c(CHRS, unique(zoom_rgns$chr[which(!zoom_rgns$chr %in% CHRS)]) ) , ordered = TRUE)
zoom_rgns = zoom_rgns[order(chr)]


#
#
#
oneprobe = probes[!duplicated(probes$probe_id)]
has_both = findOverlapPairs(toGRanges(duplicons) ,toGRanges(oneprobe)+50000 )
onedup = as.data.table(first(has_both))
onedup$probe_id = second(has_both)$probe_id
onedup$label = as.numeric(as.factor(onedup$probe_id))
probe_l = ggplot(data=onedup)+
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=.25, fill=hex)) +
  geom_segment(data= onedup %>%
                 group_by(probe_id,label) %>% 
                 summarise(start=min(start), end=max(end)),
               aes(x=start, xend=end, y=.55, yend=.55, color=probe_id), size=5, alpha=.85)+
  facet_wrap(vars(label),ncol=nrow(oneprobe), scales="free_x") + 
  scale_fill_identity()+
  scale_color_manual(values = colors)+
  theme_map()+theme(legend.position = "none",panel.border = element_rect(colour = "black", fill=NA, size=1))
probe_l
#
#
#

#colnames(CYTO.df)= c("chr","start","end","gieStain")
#x=CYTO.df
#x$end =20e6
#g = c(toGRanges(x), toGRanges(CENS))
#kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, chromosomes = CHRS[CHRS %in% c(ACHRO, "chr2","chr9")])

#pdf("acro_duplicons_flat.pdf", height = 10, width = 16)
rm.z = zoom(rm, zoom_rgns)
sedef_flat.z = zoom(sedef_flat, zoom_rgns)
duplicons.z = zoom(duplicons, zoom_rgns)
probes.z = zoom(probes, zoom_rgns)
Ns.z = zoom(NS, zoom_rgns)
#probes.o = probes.z[probes.z$origin == 1]
#x1 = toGRanges(sedef_flat.z[,c("zname","start","end")])
#x2 = toGRanges(sedef_flat.z[,c("zname2","start2","end2")])

flat_duplicon=as.ggplot(expression(
  kp<-plotKaryotype(genome=data.frame(chr=zoom_rgns$name, length=zoom_rgns$end-zoom_rgns$start)),
  #kpPlotLinks(kp, data = x1, data2 = x2,  border=NA, col=transparent("orange", 0.8), data.panel = "ideogram"),
  kpAddCytobands(kp),
  kpRect(kp, chr=duplicons.z$zname, x0 = duplicons.z$start, x1=duplicons.z$end, y0=0, y1=1, r0=0, r1=.4,  border=NA, col=duplicons.z$hex ),
  kpPlotRegions(kp, 
                data=toGRanges(probes.z[,c("zname","start","end")])+50000, 
                r0=0.5, r1=0.7, 
                col=probes.z$color, border=darker(probes.z$color, 50),
                avoid.overlapping = FALSE),
  kpPlotMarkers(kp, chr=probes.z$zname, x=(probes.z$start + probes.z$end)/2, 
                labels = probes.z$label, text.orientation = "horizontal", 
                r0=0.6, r1=.8,
                line.color = probes.z$color),
  #kpRect(kp, chr=sedef_flat.z$zname, x0 = sedef_flat.z$start, x1=sedef_flat.z$end, y0=0, y1=1,  border=NA, col=NEWCOLOR, data.panel = "ideogram"),
  kpRect(kp, chr=rm.z$zname, x0 = rm.z$start, x1=rm.z$end, y0=0.05, y1=.95,  border=NA, col=transparent(OLDCOLOR,0.25), data.panel = "ideogram"),
  kpRect(kp, chr=Ns.z$zname, x0 = Ns.z$start, x1=Ns.z$end, y0=0.05, y1=.95,  border=NA, col="black", data.panel = "ideogram")
));

probe_plot = plot_grid(probe_l, flat_duplicon, nrow=2, rel_heights = c(1,9), scale = c(0.7, 1)); probe_plot

ggsave("acro_duplicons_flat.pdf", height = length(unique(duplicons.z$zname))*.8, width = 16, plot=probe_plot)

exit()


ggplot() +
  geom_segment(data=rm.z, aes(x=start, xend=end, y=zname, yend=zname)) + 
  facet_grid(data=rm.z, cols = var(zname), rows = var(chr)) + theme_cowplot()

#dev.off()


#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")




chm13 = readbed("../Assembly_analysis/Masked/{V}.duplicons.bed", "T2T CHM13")
hg38 = readbed("data/misc_files/hg38.no_alt.duplicons.bed.gz", "GRCh38")
hg19 = readbed("data/misc_files/hg19.no_alt.duplicons.bed.gz", "GRCh37")

sedef= readbed("../Assembly_analysis/SEDEF/{V}.SDs.bed", "T2T CHM13")
low= readbed("../Assembly_analysis/SEDEF/{V}.SDs.lowid.bed", "T2T CHM13")
rm = readbed("../Assembly_analysis/Masked/{V}_repeatmasker.out.bed", "T2T CHM13", rm=T)
sat = rm[rm$type == "Satellite"]

duptaged = rgntag(chm13, sedef, "SD")
sedef = rgntag(sedef, chm13, "Duplicon")


#
#
#
#

dups = data.table(rbind(chm13, hg38) %>% 
                    group_by(duplicon,Assembly) %>% 
                    summarise(count=sum(end-start), anc= unique(ancestral_position), chr_band= unique(chr_band)) %>%
                    mutate(diff=max(count)-min(count)) %>%
                    ungroup()%>%arrange(desc(diff))
); dups

dup.genes = fread("tmp.dup.gene"); colnames(dup.genes)=c("duplicon","gene")
dups=merge(dups, dup.genes, by="duplicon", all.x=T)
dups$label = paste(dups$duplicon, paste0("(",paste0(str_trunc(dups$gene, 35),")")))
dups=data.table(dups%>%arrange(desc(diff)))
dups$label = factor(as.character(dups$label), levels = unique(as.character(dups$label)))
dups
#dups$label = paste(dups$duplicon, paste0("(", paste0(sub("NA","",dups$chr_band),")")) )



n = 30*length(unique(dups$Assembly))
duplicons = ggplot(data = dups[1:n], aes(y=label, weight=count/1000, fill=Assembly)) + 
  geom_bar(position="dodge") +
  scale_fill_manual(values=COLORS)+
  xlab("Total duplicon content (kbp)") +
  ylab("Duplicon (Gene:count,...)") +
  scale_y_discrete(position = "right")+
  scale_x_reverse(label=comma)+
  labs(title="30 duplicons with largest change between\nT2T CHM13 and GRCh38") +
  theme_cowplot()+theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5));duplicons

ggsave("duplicons.pdf", plot=plot_grid(duplicons, flat_duplicon, labels=c("a","b")), width=13*1.5, height=6*1.5) 












