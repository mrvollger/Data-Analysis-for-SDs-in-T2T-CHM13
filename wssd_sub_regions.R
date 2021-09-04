#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

maxcn=15
colors = c(`T2T CHM13`=NEWCOLOR,CHM13=NEWCOLOR, GRCh38=OLDCOLOR, Polymorphic="gray", Neither = "green", `Equal CN`="gray")

ignore = c("chr", "start","end","name","chm13", "hg38")
cols =  which(!colnames(RDG) %in% ignore)

popinfo = POPINFO 
rdg = data.table(pivot_longer(RDG, cols=cols, names_to = "Samples", values_to = "CN") %>% 
                   mutate(chm13=round(chm13), hg38=round(hg38), CN=round(CN)) %>%
                   drop_na() %>%
                   filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
                   mutate(chm13_diff = (CN-chm13), hg38_diff=(CN-hg38)) %>% 
                   mutate(diff = case_when(abs(chm13_diff) <= abs(hg38_diff) ~ chm13_diff, abs(chm13_diff) > abs(hg38_diff) ~ hg38_diff), 
                          winner = case_when(abs(chm13_diff)==abs(hg38_diff) ~ "Equal CN", abs(chm13_diff) < abs(hg38_diff) ~ "CHM13",  abs(chm13_diff) > abs(hg38_diff) ~ "GRCh38"))
                   )


# remove entries that overlap with Ns 
has_Ns = queryHits(findOverlaps(toGRanges(rdg[,c("chr","start","end","name")]), toGRanges(NS)))
rdg = rdg[-has_Ns,]

rdg$color = colors[rdg$winner]
n.sms = length(unique(rdg$Samples))


# make a region df that defines winner for each region
min_n_to_support = 0
n_stds = 2
rgn.df  = data.table(rdg %>% group_by(chr,start,end,name) %>%
                        summarise( CHM13 = sum( abs(CN-chm13) <= 1 ), GRCh38=sum(abs(CN-hg38) <= 1), 
                                   mean=mean(CN), median=median(CN), std = sd(CN),
                                   chm13_cn =unique(chm13)[1], hg38_cn = unique(hg38)[1]) %>%
                        mutate( status = case_when( 
                                                    CHM13==0 & GRCh38==0 ~ "Neither", 
                                                    CHM13 >  GRCh38  ~ "CHM13", 
                                                    GRCh38 > CHM13 ~ "GRCh38",
                                                    TRUE ~ "Equal CN"),
                                region_color = colors[status],
                                status_std = case_when( 
                                  abs(chm13_cn - median) > std * n_stds & abs(hg38_cn - median) > std * n_stds ~ "Neither", 
                                  abs(chm13_cn - median) <= std * n_stds & abs(hg38_cn - median) <= std *n_stds ~ "Polymorphic", 
                                  abs(chm13_cn - median) <= std * n_stds ~ "CHM13", 
                                  abs(hg38_cn - median) <= std * n_stds ~ "GRCh38"
                                  ),
                                region_color_std = colors[status_std]
                                )
                     )
                        

# remove very high CN regions to remove the rDNA decoy copies
to_remove_for_high_cn_and_rdna = rgn.df[median > 200 & chr %in% ACHRO]$name
rgn.df = rgn.df[! name %in% to_remove_for_high_cn_and_rdna]
rdg = rdg[! name %in% to_remove_for_high_cn_and_rdna]




# merge wiht super pop info and location
pop_location = data.table(super_pop = c("AMR","AFR","EA","SA","SIB","WEA","OCN"),
                          lat  = c( 13, 7.7, 33, 10.3, 61, 52, 1),
                          long  = c(-88, 21,  110, 75.5, 95, 27, 130)
); 
pop_location$lat = as.numeric(pop_location$lat); pop_location$long = as.numeric(pop_location$long)
rdg$sample = rdg$Samples
rdg=merge(rdg, popinfo, by="sample")
rdg= data.table( rdg %>% merge(pop_location, by="super_pop") %>% group_by(super_pop)%>% mutate(n_super_pop = length(unique(sample)))) 




# importatn TBC1D3
#ggplot(data=rdg %>% filter(name == "chr17_38963566_38978236"))+geom_count(aes(x=CN, y=1))



# merge with region info
rdg = merge(rdg, rgn.df[,c("name","status")], on="name")
rdg$region_color = colors[rdg$status]

#
#
# This plot shows the cumlative distribution of CN imporvment
#
rdg_perfect = copy(rdg)
rdg_perfect = merge(rdg, rgn.df[,c("name","mean", "median")], on="name")
rdg_perfect$perfect = rdg_perfect$median
rdg_perfect$perfect_diff = rdg_perfect$CN - rdg_perfect$median

ecdf_rdg = pivot_longer(rdg_perfect, cols = which(colnames(rdg_perfect) %in% c("chm13_diff","hg38_diff", "perfect_diff")),
                        names_to="ecdf_name", values_to = "ecdf_diff" )

simple_auc = function(fecdf, maxval){
  auc = sum(fecdf(seq(0, maxval)))/maxval
  auc
}
vals13 = abs( ecdf_rdg[ecdf_rdg$ecdf_name =="chm13_diff",]$ecdf_diff)
ecdf13 = ecdf(vals13)
vals38 = abs( ecdf_rdg[ecdf_rdg$ecdf_name =="hg38_diff",]$ecdf_diff)
ecdf38 = ecdf(vals38)
vals_per = abs( ecdf_rdg[ecdf_rdg$ecdf_name =="perfect_diff",]$ecdf_diff)
ecdf_per = ecdf(vals_per)

max_allowed_cn_diff = 30
auc.df = data.frame(Reference = c("T2T CHM13", "GRCh38", "Median CN of SGDP"), 
                    AUC = round(c(simple_auc(ecdf13, max_allowed_cn_diff), 
                            simple_auc(ecdf38, max_allowed_cn_diff),
                            simple_auc(ecdf_per, max_allowed_cn_diff)
                            ),
                            2)
                    )
#                    color = c("chm13_diff", "hg38_diff"))
auc.df
ecdf_rdg$alpha = 1; ecdf_rdg[ecdf_rdg$ecdf_name == "perfect_diff",]$alpha = 0.9

ecdf_cn_plot = ggplot(data=ecdf_rdg) + 
  stat_ecdf(data = . %>% filter(ecdf_name != "perfect_diff"),
            aes(x=abs(ecdf_diff), color=ecdf_name), size=2, alpha=.9) +
  stat_ecdf(data = . %>% filter(ecdf_name == "perfect_diff"),
            aes(x=abs(ecdf_diff), color=ecdf_name), alpha = 0.75,  size=1) +
  scale_color_manual(values=c(chm13_diff=NEWCOLOR, hg38_diff=OLDCOLOR, perfect_diff="forestgreen"),
                     labels=c("T2T CHM13", "GRCh38", "Median CN of SGDP"))+
  #coord_cartesian(xlim = c(0,20))+
  facet_zoom(x= abs(ecdf_diff) <= max_allowed_cn_diff, zoom.size = 3)+
  xlab("Allowed CN difference between\nthe sample and the reference")+
  ylab("Cumulative fraction of correct CN representation")+
  theme_cowplot()+theme(legend.position = "top", 
                        legend.title = element_blank())

ecdf_cn_plot = ggdraw() + 
  draw_plot(ecdf_cn_plot) +
  draw_grob(tableGrob(auc.df, theme=ttheme_minimal(), rows=NULL), 
            x=1, y=1, hjust = 1, vjust = 1.25, halign = 1, valign = 1)
ecdf_cn_plot

ggsave(glue("{SUPP}/wssd_ecdf.pdf"), width = 12, height = 8, plot=ecdf_cn_plot)


#
# This figure shows the CN of chm13 and hg38 vs the median for each loci
#

#ggplot(data = rgn.df %>%
#         arrange(median) ) + 
#  geom_link(aes(x=name, xend=name, y=chm13_cn, yend=median))



#
#
#
p = ggplot(data=rdg%>% mutate(diff=case_when(diff>maxcn ~ maxcn, diff < -maxcn ~ -maxcn, TRUE ~ diff)))+
  geom_histogram(aes(x=diff, fill=winner), binwidth = 1, position=position_dodge(0.8)) + 
  xlab("Copy number difference (sample CN - reference CN)") + ylab("# of sample genotypes with observed CN difference") +
  scale_fill_manual(values=colors) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels =  c(glue("-{maxcn}+"),(1-maxcn):(maxcn-1), glue("{maxcn}+") ), breaks = -maxcn:maxcn) +
  coord_cartesian(xlim=c(-maxcn,maxcn))+ theme_cowplot() + theme(legend.position = "none");p
ggsave(glue("{SUPP}/wssd_cn_diff_hist.pdf"), width = 12, height = 8, plot=p)


#
# wssd karyoplot
#
pp=plotDefaultPlotParams(plot.type = 2); 
pp$plot.params$ideogramheight = 200
ideo = as.ggplot(expression(
  kp <- plotKaryotype(genome=GENOME, chromosomes = NOYM[NOYM %in% unique(rdg$chr)], cytobands = CYTO, plot.type = 2, plot.params = pp$plot.params),
  topp<-rdg[rdg$winner %in% c("CHM13") & abs(rdg$chm13_diff) <= 1 ,],
  botp<-rdg[rdg$winner %in% c("GRCh38") & abs(rdg$hg38_diff) <= 1,],
#kpPlotDensity(kp, data = toGRanges(rdg[winner=="Polymorphic",1:3]), col="gray")
kpPlotRegions(kp, data = toGRanges(topp[,c("chr","start","end")]), col=topp$color, ylim=n.sms, data.panel = 1, r0=.15),
kpPlotRegions(kp, data = toGRanges(botp[, c("chr","start","end")]), col=botp$color, ylim =n.sms, data.panel = 2, r0=.15),
kpPlotRegions(kp, 
              data = toGRanges(rgn.df[rgn.df$status!="Equal CN",c("chr","start","end")]), 
              col=rgn.df$region_color[rgn.df$status!="Equal CN"], 
              ylim = n.sms, data.panel = "ideogram", r0=0.2, r1=0.8, border = NA)
)); ideo
ggsave(glue("{SUPP}/wssd_ideo_exact_hist.pdf"), width = 12, height = 8, plot=ideo)


#
# winner plot
#
t1 = rdg %>% group_by(winner) %>% summarise(count=length(winner))
t2 = rgn.df %>% group_by(status) %>% summarise(count = length(status))

p3 = ggplot(data=t1) + geom_bar(aes(x=winner, y=count, fill=winner), stat = "identity") + 
  geom_text(aes(label=comma(count), x=winner, y=count), vjust=-1, angle=0) +
  scale_fill_manual(values=colors) + theme_cowplot() + 
  scale_y_continuous(label=comma)+
  ylab("# of sample genotypes with CN closer to CHM13 or GRCh38") + xlab("") +
  theme(legend.position = "none") ;p3
ggsave(glue("{SUPP}/wssd_sample_closer_hist.pdf"), width = 8, height = 8, plot=p3)

p4 = ggplot(data=t2) + geom_bar(aes(y=status, x=count, fill=status), stat = "identity") + 
  geom_text(aes(label=count, y=status, x=count), vjust=-1, angle=0) +
  scale_fill_manual(values=colors) + theme_cowplot() + xlab(glue("# of regions where more samples genotype exactly with either CHM13 or GRCh38")) + ylab("") + theme(legend.position = "none"); p4




#
# world map
#
sizescale= 2
zz_extra = data.table(rdg %>% group_by(super_pop,lat,long,n_super_pop,chr,start,end,name) %>%
                        summarise( CHM13 = sum( abs(CN-chm13) <= 1 ), GRCh38=sum(abs(CN-hg38) <= 1)) %>%
                        mutate( status = case_when( CHM13==0 & GRCh38==0 ~ "Neither", CHM13 >  GRCh38  ~ "CHM13", GRCh38 > CHM13 ~ "GRCh38", TRUE ~ "Equal CN" )) %>%
                        ungroup() %>% group_by(super_pop,lat,long,n_super_pop) %>%
                        summarise(CHM13 = sum(status=="CHM13")/length(status), GRCh38 = sum(status=="GRCh38")/length(status),`Equal CN` = sum(status=="Equal CN")/length(status), radius=sizescale*sqrt(unique(n_super_pop))));#unique(n_super_pop)/sizescale));zz_extra



world <- map_data('world')
map <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill=NA, color="black") +
  coord_quickmap() +
  geom_scatterpie(aes(x=long, y=lat, group=super_pop, r=radius),
                     data=zz_extra, cols=c("CHM13", "GRCh38", "Equal CN"), color=NA, alpha=0.95) +
  geom_scatterpie_legend(zz_extra$radius, x=-160, y=-55, labeller = function(x){x*x/sizescale} ) +
  scale_fill_manual(values=colors) + theme_map() + ylab("Latitude") + xlab("longitude") + theme(legend.position = "top", legend.title = element_blank(), legend.justification = "center"); map

ggsave(glue("{SUPP}/wssd_map.pdf"), width = 12, height = 8, plot=map)


#
# mbp histogram
# 
bp_hist = ggplot(data=rgn.df,
                 aes(x=status, weight=(end-start)/1e6, fill=status)
                 )+
  geom_bar()+
  theme_cowplot()+
  scale_y_continuous(labels = comma) + 
  scale_fill_manual(values=colors) +
  geom_text(stat='count', aes(label=round(..count..,2)), vjust=-1)+
  xlab("") + ylab("# of Mbp where more samples have CN\ngenotypes equal to the reference") +
  theme(legend.position = "none")
bp_hist
ggsave(glue("{SUPP}/wssd_bp_hist.pdf"), width = 8, height = 8, plot=bp_hist)

bp_hist_std = ggplot(data=rgn.df,
       aes(x=status_std, weight=end-start, fill=status_std) 
       )+
  geom_bar()+
  theme_cowplot()+
  scale_y_continuous(labels = comma) + 
  scale_fill_manual(values=colors) +
  geom_text(stat='count', aes(label=comma(..count..)), vjust=-1)+
  xlab("") + ylab(glue("# of bp where the reference is within the median CN +/- {n_stds} sd")) +
  theme(legend.position = "none")
bp_hist_std
ggsave(glue("{SUPP}/wssd_bp_hist_std.pdf"), width = 8, height = 8, plot=bp_hist_std)

