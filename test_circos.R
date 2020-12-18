#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")
sedef= SEDEF
duplicons = DM_BED #fread("../Assembly_analysis/Masked/chm13.draft_v1.0_plus38Y.duplicons.bed")#DM
enriched = readbed("../Assembly_analysis/SEDEF/{V}.sedef.enriched.bed", "highsd")[,1:4]
cyto.df = data.frame(CYTO)[, c("seqnames","start","end","seqnames","gieStain")]
set.seed(3)
chr_bg_color = rand_color(length(cyto.df$seqnames), transparency = 0.6)
names(chr_bg_color) =cyto.df$seqnames

plotname="circos.pdf"

if("acro"=="acrox"){
  #cyto.df = cyto.df[cyto.df$seqnames %in% ACHRO,]
  #cyto.df$end[cyto.df$end > 20e6]=20e6
  sedef=sedef[chr %in% ACHRO | chr2 %in% ACHRO]
  highest_achro = sedef %>% 
    mutate(xx=pmin(as.character(chr), as.character(chr2)), yy = pmax(as.character(chr), as.character(chr2))) %>%
    group_by(xx,yy) %>% 
    summarise(total=sum(alnB)) %>% 
    arrange(-total) %>% 
    head(20)
    #filter(total>1e6)
  names=unique(c(highest_achro$xx, highest_achro$yy));names;length(names)

  sedef=sedef[chr %in% names & chr2 %in% names]
  cyto.df=cyto.df[cyto.df$seqnames %in%  names,]
  duplicons=duplicons[chr %in% names]
  enriched = enriched[chr %in% names & end-start > 100000]
  plotname="circos_acro_enriched.pdf"
}else if ("achrox2" =="achro2"){
  sedef=SEDEF[chr %in% ACHRO | chr2 %in% ACHRO]
  highest_achro = sedef %>% 
    mutate(xx=pmin(as.character(chr), as.character(chr2)), yy = pmax(as.character(chr), as.character(chr2))) %>%
    group_by(xx,yy) %>% 
    summarise(total=sum(alnB)) %>% 
    arrange(-total) %>% 
    head(30)
  names=unique(c(highest_achro$xx, highest_achro$yy));names;length(names)
  #names=c("chr9", "chr13")
  sedef = sedef[chr %in% names & chr2 %in% names & !(chr %in% ACHRO & chr2 %in% ACHRO) & (chr %in% ACHRO | chr2 %in% ACHRO) & alnB>1000 & fracMatch>0.8 & chr!=chr2]
  duplicons=duplicons[chr %in% names]
  does_overlap = overlaps(duplicons, sedef)
  duplicons_overlap = duplicons[does_overlap]
  duplicons$overlap = does_overlap
  
  #duplicons$hex = sapply(strsplit(duplicons$V9, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
  
  enriched = enriched[chr %in% names & end-start > 1000000]
  enriched = readbed("large.enriched.bed", "x")[,1:4]
  enriched_sds = rgntag(sedef, enriched, "enriched"); enriched_sds=enriched_sds[enriched==TRUE]
  
  zoom_rgns  = data.table(enriched_sds %>% 
                            group_by(chr) %>% 
                            summarise(start=min(start), end=max(end)) %>% mutate(V4=(end-start))) #, label=paste0(chr,":",start,"-",end)) )
  enriched = enriched[chr %in% names]
  #enriched = zoom_rgns
  sedef=sedef[chr %in% names & chr2 %in% names]
  cyto.df=cyto.df[cyto.df$seqnames %in%  names,]
  duplicons=duplicons[chr %in% names]
  enriched = enriched[chr %in% names & end-start > 100000]
  plotname="circos_acro_enriched.pdf"
  dim(duplicons)
  dim(sedef)
  dim(enriched)
  dim(duplicons_overlap)
  enriched
} else if("fast"=="faxst"){
  N=1000
  duplicons=duplicons[sample(.N,N)]
  sedef=sedef[sample(.N,N)]
  plotname="circos_test.pdf"
} else{
  enriched=enriched[end-start > 1000000]
}

cyto.df$seqnames =as.character(cyto.df$seqnames)
cyto.df




#sedef=sedef[paste(sedef$chr,sedef$start1) >= paste(sedef$chr2, sedef$start2) & sedef$aln_len > 10000 ]
#sedef=sedef[sedef$aln_len > 25000 ]
sedef=sedef[order(aln_len)]
print(dim(sedef))
print(sum(paste(sedef$chr,sedef$start1) > paste(sedef$chr2, sedef$start2)))

#sedef=data.table(sample_n(sedef,1000))

sedef.inter=sedef[chr != chr2 ]
sedef.intra=sedef[chr == chr2 & chr > chr2]


new = overlap_either(sedef.inter, NEW)
sedef.inter$color = GRAY
sedef.inter$color[new] = NEWCOLOR


nrow(sedef.inter)
nrow(sedef.intra)
#sedef=achro(sedef)
#sedef=sample_n(sedef[chr!=chr2], 1000)

sdcolor=function(df, amount=0.5){
  color=rep("#000000", nrow(df))
  color[df$aln_len >= 10000 & df$fracMatch > 0.95]="#FF0000"
  return(transparent(color, amount = amount))
}

convert_data = function(bed, conversion){
  which_region <- function(row) {
    matches = conversion[row["chr"] == conversion$chr & row["end"] >= conversion$start & row["start"] <= conversion$end ]; 
    label=matches$label[1]
    return(label)
  }
  #bed$label = apply(bed[, 1:3], 1, which_region)
  #bed = bed[!is.na(bed$label)]
  tmp = data.table(data.frame(second(do_bedtools_intersect(toGRanges(bed), toGRanges(conversion[,4:6],f=1.0),  loj = TRUE))))
  cond= tmp$end != -1
  keep = tmp[cond]
  bed=bed[cond]
  bed$label=paste(paste(keep$seqnames, keep$start, sep=":"), keep$end, sep="-")
  return(bed)
}

f1 = function(){
  gap.degree=360/(4*length(unique(cyto.df$seqnames))); gap.degree
  circos.par(cell.padding = c(0, 0, 0, 0), gap.degree=gap.degree)
  circos.initializeWithIdeogram(cyto.df, ideogram.height = .05 )
  #circos.genomicRect(enriched[,1:3], ybottom =0, ytop=1, border = NA, col="red")
  #circos.genomicRainfall(sedef[,1:3])
  #circos.genomicRainfall(sedef.inter[,1:3]);  circos.genomicRainfall(sedef.inter[,c("chr2","start2", "end2")])
  circos.genomicLink(sedef.inter[,c("chr","start", "end")], sedef.inter[,c("chr2","start2", "end2")], border=NA, col=sedef.inter$color)#col=sdcolor(sedef.inter, amount = 0.35))
};
f2 = function(){
  gap.degree=360/(4*length(conversion$label)); gap.degree
  circos.par(cell.padding = c(0, 0, 0, 0), track.height=0.25/2, gap.degree=gap.degree)
  circos.initialize(factors=conversion$label, xlim = conversion[,c(2,3)])
  color = sapply(strsplit(zoomed_data$V9, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255) )
  color_2 = sapply(strsplit(zoomed_data_2$V9, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255) )
  
  YMAX = max(zoomed_data$end)
  # plot the track.
  circos.track(zoomed_data$label, x = zoomed_data$start, y = zoomed_data$end,
              panel.fun = function(x, y) {
                l = zoomed_data$label == CELL_META$sector.index
                circos.rect(x, 2, y, YMAX, col=color[l], border=NA)
              }
  )
  # plot the sub track
  circos.track(zoomed_data_2$label, x = zoomed_data_2$start, y = zoomed_data_2$end,
               panel.fun = function(x, y) {
                 l = zoomed_data_2$label == CELL_META$sector.index
                 circos.rect(x, 2, y, YMAX, col=color_2[l], border=NA)
               }
  )
}

enriched$label=paste(paste(enriched$chr, enriched$start, sep=":"), enriched$end, sep="-")
conversion=unique(enriched[,c("label","start","end","chr")])
conversion$label=as.character(conversion$label)
conversion$start2=conversion$start
conversion$end2 = conversion$end
zoomed_data = convert_data(duplicons, conversion)
zoomed_data_2 = convert_data(duplicons_overlap, conversion)

DIM = 7
pdf(paste0("simple_", plotname), height = DIM, width =DIM)
f1()
dev.off()

pdf(plotname, height =DIM, width = DIM)
circos.nested(f2, f1, conversion, connection_col =chr_bg_color[conversion[[4]]], connection_height=0.3, connection_border = NA)
dev.off()
dev.off()



exit

sedef.inter = readbed("../Assembly_analysis/SEDEF/tmp.bed", "tmp")[chr != chr2]; sedef.inter
sedef.inter$color = NEWCOLOR
f1()

pairs = data.table(read_excel("../Assembly_analysis/SEDEF/{V}.sedef.summary.xlsx", sheet="bp"))
pairs38 = data.table(read_excel("../Assembly_analysis/SEDEF/hg38.chr_only.sedef.summary.xlsx", sheet="bp"))

x = as.matrix(pairs[,-1])
row.names(x) = pairs$...1
heat_df = melt(x)

x = as.matrix(pairs38[,-1])
row.names(x) = pairs38$...1
heat38_df = melt(x)

heat_df$hg38 = heat38_df$value
heat_df$diff = heat_df$value - heat_df$hg38
heat_df$log2 =log2(heat_df$value / heat_df$hg38)

heat_df$label = round(heat_df$value/10^6,2)
#heat_df$label[heat_df$label < 1 ]=NA
heat_df$cut = cut(heat_df$label, breaks=c(quantile(heat_df$label, probs = seq(0, 1, by = 1/9))), include.lowest=TRUE)

p = ggplot(data=heat_df, aes(x=Var1, y=Var2, fill=cut, label=label)) +
  geom_tile() + 
  geom_text(color="white") +
  theme_cowplot()+ xlab("Destination chromosome")+ylab("Source chromosome")+ggtitle("Non-redundant Mbp on destination chromosome")+
  theme(legend.position = "bottom", axis.text.x = element_text(angle=45, hjust = 1))+ coord_fixed()+
  scale_fill_brewer(palette = "YlOrRd", direction = 1);p

p2=ggplot(data=heat_df, aes(x=value, fill=cut, label=label)) +
  geom_histogram(bins=20, color="black")+ 
  scale_x_continuous(trans="log10", labels = comma)+
  theme_cowplot()+ xlab("bp of SD")+
  theme(legend.position = "none")+
  scale_fill_brewer(palette = "YlOrRd", direction = 1)

heat_df$label = round(heat_df$diff/10^6,2)
#heat_df$label[heat_df$label < 1 ]=NA
heat_df$cut = cut(heat_df$label, breaks=c(quantile(heat_df$label, probs = seq(0, 1, by = 1/9))), include.lowest=TRUE)
p3= ggplot(data=heat_df, aes(x=Var1, y=Var2, fill=log2, label=log2)) +
  geom_tile() + 
  geom_text(color="white") +
  theme_cowplot()+ xlab("Destination chromosome")+ylab("Source chromosome")+ggtitle("Non-redundant Mbp on destination chromosome");p3#+
  theme(legend.position = "bottom", axis.text.x = element_text(angle=45, hjust = 1))+ coord_fixed()+
  scale_fill_brewer(palette = "RdBu", direction = -1);p3

ggsave("SD_heat_map.pdf", plot=plot_grid(p2,p, nrow=2, rel_heights = c(1,4)), height = 16, width = 12)

