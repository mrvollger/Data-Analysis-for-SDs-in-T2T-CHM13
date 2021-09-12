odir = "~/mvollger@uw.edu - Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/added-after-2021-09-04/structured-abstract/"
sds13 = SEDEF_V1.1[SEDEF_V1.1$original & (acro | acro2) & aln_len > 1e4]
sds38 = SEDEF_38[SEDEF_38$original  & (acro | acro2) & aln_len > 1e4 ]


dim(sds13)
dim(sds38)
chrs_for_cmp = unique(c(
  as.character(sds13$chr),
  as.character(sds13$chr2),
  as.character(sds38$chr),
  as.character(sds38$chr2)
  ))
tmp = sds13 %>% group_by(chr,chr2) %>% summarise(n=n()) %>% ungroup() %>% group_by(chr) %>% summarise(n = sum(n)) %>%
  arrange(-n) %>% head(10)
tmp
chrs_for_cmp = unique(as.character(tmp$chr))
length(chrs_for_cmp)

#
#
#
chr_bg_color = rand_color(length(chrs_for_cmp), transparency = 0.6)
#
#
#
enriched = ENRICHED[chr %in% chrs_for_cmp & overlaps(ENRICHED, sds13)] 
# merge them
e = reduce(toGRanges(enriched))
g = gaps(e)
enriched = as.data.table(reduce(c(e,g[width(g)< 1e6])))
colnames(enriched)[1]="chr"
dim(enriched)

duplicons = DM_BED[overlaps(DM_BED,enriched)]
duplicons_38 =  DM_BED_38[overlaps(DM_BED_38, enriched)]
N_38 = as.data.table(intersect(toGRanges(Ns_38), toGRanges(enriched)))
N_38$V9 = "0,0,0"; colnames(N_38)[1] = "chr"
N_38$chr = factor(N_38$chr, levels = CHRS, ordered = TRUE)
duplicons_38 = bind_rows(duplicons_38, N_38[,c("chr","start","end","V9")]) %>% data.table()
duplicons_38 =  duplicons_38[overlaps(duplicons_38, enriched)]

convert_data = function(bed, conversion){
  which_region <- function(row) {
    matches = conversion[row["chr"] == conversion$chr & row["end"] >= conversion$start & row["start"] <= conversion$end ]; 
    label=matches$label[1]
    return(label)
  }
  #bed$label = apply(bed[, 1:3], 1, which_region)
  #bed = bed[!is.na(bed$label)]
  tmp = data.table(data.frame(
    second(do_bedtools_intersect(toGRanges(bed), toGRanges(conversion[,4:6], f=1.0),  loj = TRUE))
    ))
  cond= tmp$end != -1
  keep = tmp[cond]
  bed=bed[cond]
  bed$label=paste(paste(keep$seqnames, keep$start, sep=":"), keep$end, sep="-")
  return(bed)
}
enriched$label=paste(paste(enriched$chr, enriched$start, sep=":"), enriched$end, sep="-")
conversion=unique(enriched[,c("label","start","end","chr")])
conversion$label=as.character(conversion$label)
conversion$start2=conversion$start
conversion$end2 = conversion$end
zoomed_data_13 = convert_data(duplicons, conversion)# %>% sample_n(1000) %>% data.table()
zoomed_data_38 = convert_data(duplicons_38, conversion)# %>% sample_n(1000) %>% data.table()

dim(zoomed_data_13)
dim(zoomed_data_38)
#
#
#


circle <- function(color=NEWCOLOR){
  # data
  df.a = df.a[chr %in% chrs_for_cmp & chr2 %in% chrs_for_cmp]
  sedef.inter <- df.a[chr != chr2]# %>% sample_n(500) %>% data.table()
  sedef.intra <- df.a[chr == chr2] 
  print(nrow(sedef.inter))
  
    
  cyto.df <- data.frame(CYTOFULL)[, c("seqnames","start","end","seqnames","gieStain")]
  cyto.df = cyto.df[cyto.df$seqnames %in% chrs_for_cmp,]
  #cyto.df$seqnames = factor(cyto.df$seqnames, levels =  c(CHRS[!CHRS %in% ACHRO ], ACHRO))
  cyto.df  = cyto.df[order(cyto.df$gieStain, cyto.df$seqnames),]
  cyto.df$seqnames = as.character(cyto.df$seqnames)
  
  #new <- overlap_either(sedef.inter, NEW)
  #sedef.inter$color = GRAY
  #sedef.inter$color[new] = NEWCOLOR
  sedef.inter$color = color
  # ploting
  gap.degree=360/(8*length(unique(cyto.df$seqnames))); gap.degree

  circos.par(cell.padding = c(0, 0, 0, 0), 
             gap.degree=gap.degree
             )
  circos.initializeWithIdeogram(cyto.df,ideogram.height = .10)
  #circos.genomicRect(enriched[,1:3], ybottom =0, ytop=1, border = NA, col="red")
  #circos.genomicRainfall(sedef[,1:3])
  #circos.genomicRainfall(sedef.inter[,1:3]);  circos.genomicRainfall(sedef.inter[,c("chr2","start2", "end2")])
  circos.genomicLink(sedef.inter[,c("chr","start", "end")],
                     sedef.inter[,c("chr2","start2", "end2")],
                     border=NA, col=sedef.inter$color)
  
  print(nrow(sedef.intra))
  circos.genomicLink(sedef.intra[,c("chr","start", "end")],
                     sedef.intra[,c("chr2","start2", "end2")],
                     border=NA,
                     col=NEWCOLOR
                     )
}

f2 = function(){
  gap.degree=360/(6*length(conversion$label)); gap.degree
  length(unique(conversion$label))
  circos.par(cell.padding = c(0, 0, 0, 0),
             gap.degree = gap.degree,
             #"start.degree" = -55,
             track.height=0.25/4)
  
  circos.initialize(factors=conversion$label, xlim = conversion[,c(2,3)])
  
  color = sapply(strsplit(zoomed_data$V9, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255) )
  
  YMAX = max(zoomed_data$end)
  # plot the track.
  circos.track(zoomed_data$label, x = zoomed_data$start, y = zoomed_data$end,
               panel.fun = function(x, y) {
                 l = zoomed_data$label == CELL_META$sector.index
                 circos.rect(x, 2, y, YMAX, col=color[l], border=NA)
               }
  )
}


pdf(glue("{odir}/grch38.pdf"), height = 6, width = 6)
df.a <<- sds38
zoomed_data <<- zoomed_data_38
circos.nested(f2, circle, conversion, connection_col=chr_bg_color[conversion[[4]]], connection_height=0.1, connection_border = NA)
dev.off()

pdf(glue("{odir}/t2t.pdf"), height = 6, width = 6)
df.a <<- sds13
zoomed_data <<- zoomed_data_13
circos.nested(f2, circle, conversion, connection_col=chr_bg_color[conversion[[4]]], connection_height=0.1, connection_border = NA)
dev.off()

pdf(glue("{odir}/grch38_no_dup.pdf"), height = 6, width = 6)
df.a <<- sds38
circle()
dev.off()

pdf(glue("{odir}/t2t_no_dup.pdf"), height = 6, width = 6)
df.a <<- sds13
circle()
dev.off()


####################################################################################################
####################################################################################################
####################################################################################################


offlong = -24

# read pop info
sgdp_pop_info = fread("http://simonsfoundation.s3.amazonaws.com/share/SCDA/datasets/10_24_2014_SGDP_metainformation_update.txt",
                      blank.lines.skip=T,
                      drop=seq(12,24)
) %>% drop_na() %>%
  group_by(Population_ID) %>%
  summarise(lat = unique(Latitude),
            long = unique(Longitude),
            Region = unique(Region)
  ) %>%
  mutate(pop = Population_ID) %>% 
  bind_rows(data.frame(
    pop = c( "Norwegian","Saami","Crete"),
    lat = c(52, 13, 52),
    long = c(27, -88, 27),
    Region = c("WestEurasia","Africa", "WestEurasia")
  )) %>%
  ungroup() %>%
  group_by(Region) %>%
  mutate(
    long2 = ifelse(long < offlong, long + 360, long), 
    long.super = mean(long2),
    lat.super = mean(lat)
    ) %>%
  dplyr::select(-Population_ID) %>% 
  distinct(pop, .keep_all = TRUE) %>%
  data.table(); sgdp_pop_info

if(F){
supers <- sort(unique(sgdp_pop_info$Region))
names(supers) = c("AFR", "AMR", "SIB", "EA", "OCN", "SA", "WEA")
sgdp_pop_info$super_pop <- names(supers)[match(sgdp_pop_info$Region, supers)]
sgdp_pop_info=sgdp_pop_info %>% 
  merge(pop_location, by="super_pop", suffixes=c("", ".super"))
}
#
# load first 70 lines for wssd_sub_regions
#

pop_places = rdg %>%
  dplyr::select(-lat, -long) %>%
  merge(sgdp_pop_info, by = "pop") %>%
  distinct(Samples, .keep_all = TRUE) %>% 
  group_by(pop, lat, long, lat.super, long.super, long2) %>%
  summarise(n_pop = n()) %>%
  mutate(xend = 200,
         yend = 10)

summary_pie_df = rdg %>%
  group_by(chr, start, end, name) %>%
  summarise(
    CHM13 = sum(abs(CN - chm13) <= 1),
    GRCh38 = sum(abs(CN - hg38) <= 1)
  ) %>%
  mutate(
    status = case_when(
      CHM13 == 0 & GRCh38 == 0 ~ "Neither",
      CHM13 > GRCh38 ~ "CHM13",
      GRCh38 > CHM13 ~ "GRCh38",
      TRUE ~ "Equal CN"
    ),
    weight = end - start
  ) %>%
  ungroup() %>%
  summarise(
    `T2T-CHM13` = sum((status == "CHM13") * weight) ,
    GRCh38 = sum((status == "GRCh38") * weight) , 
    `Equal CN` = sum((status == "Equal CN") * weight) , 
    radius = 35
  ) %>%
  mutate(lat = unique(pop_places$yend),
         long = unique(pop_places$xend)
         ) %>%
  data.table()


colors=c(colors, `T2T-CHM13`=NEWCOLOR)


#pop_places$long.super2 = ifelse(pop_places$long.super < offlong, pop_places$long.super + 360, pop_places$long.super) 
summary_pie_df$long2 = ifelse(summary_pie_df$long < offlong, summary_pie_df$long + 360, summary_pie_df$long) 

super_pop_places = pop_places %>%
  group_by(long.super, lat.super, xend, yend) %>%
  summarise(n_pop = sum(n_pop))


world <- map_data('world', wrap=c(offlong, 360 + offlong), ylim=c(-55,75))
world$long2 = world$long

map= ggplot(world, aes(long2, lat)) +
  geom_map(map=world, aes(map_id=region), fill="gray60", color="gray60") +
  #mapWorld+
  coord_quickmap(ylim = c(-50,NA)) +
  geom_curve(data = pop_places,
    aes(xend=jitter(long.super, amount=0.001), yend = lat.super, group=pop),
    alpha= 0.5,
    curvature = -0.1,
    #lineend="round",
    linetype=8,
    size = 0.5,
    color= "black"
  ) +
  geom_point(data=pop_places, 
                aes(group=pop),
                size = 1,
                color="black", fill="black"
             ) +
  geom_curve(data = super_pop_places,
             aes(
               x=long.super,
               xend=jitter(xend, amount=0.001),
               y = lat.super,
               yend=yend,
               size=log2(n_pop)
               ),
             alpha= 0.8,
             curvature = -0.7,
             arrow = arrow(angle = 5, length = unit(0.15, "inches"), type = "closed"),
             lineend="round",
             color= "black",
             angle = 130,
             ncp = 20
  ) +
  scale_size(range=c(1, 2)) +
  geom_scatterpie(
    data = summary_pie_df,
    aes(x=long2, y=lat, group=pop, r=radius),
    alpha=0.991,
    cols=c("T2T-CHM13", "GRCh38"),
    color="gray40"
  ) +
  scale_fill_manual(values=colors) +
  theme_map() +
  ylab("Latitude") + xlab("longitude") +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.justification = "center",
        plot.title = element_text(hjust=0.5)); map


ggsave(glue("{odir}/world_map_simple.pdf"), height = 6, width = 12, plot=map)







