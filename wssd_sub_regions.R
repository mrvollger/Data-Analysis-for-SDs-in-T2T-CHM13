#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")


colors = c(Polymorphic="gray", CHM13=NEWCOLOR, GRCh38=OLDCOLOR, Polymorphic="gray", Neither = "green")

tmp1 = fread("../wssd/v1.0/combined_tables/wssd/SYNTENY_EXP_SD_wssd_GMM.bed")
#tmp2 = fread("../wssd/v1.0/rdg_20kb/HIGH_CN/HIGH_CN.hgdp.combined.wssd.GMM.bed")
tmp2 = fread("../wssd/v1.0/combined_tables/wssd/HIGH_CN_wssd_GMM.bed")
if(F){
  tmp=tmp1
  maxcn = 5
  pname="low"
}else if(T){
  tmp=rbind(tmp1,tmp2) # do not do, regions tend to overlap
  maxcn=10
  pname="both"
}else{
  tmp=tmp2
  maxcn=10
  pname="high"
}

ignore = c("chr", "start","end","name","chm13", "hg38")
cols =  which(!colnames(tmp) %in% ignore) #colnames(tmp)[!colnames(tmp) %in% ignore]

popinfo = fread("../wssd/v1.0/rdg_20kb/hgdp_manifest.txt")
rdg = data.table(pivot_longer(tmp, cols=cols, names_to = "Samples", values_to = "CN") %>% 
                   drop_na() %>%
                   filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
                   mutate(chm13_diff = (CN-chm13), hg38_diff=(CN-hg38)) %>% 
                   mutate(diff = case_when(abs(chm13_diff) <= abs(hg38_diff) ~ chm13_diff, abs(chm13_diff) > abs(hg38_diff) ~ hg38_diff), 
                          winner = case_when(abs(chm13_diff)==abs(hg38_diff) ~ "Polymorphic", abs(chm13_diff) < abs(hg38_diff) ~ "CHM13",  abs(chm13_diff) > abs(hg38_diff) ~ "GRCh38"))
                   )
dim(rdg)
rdg$color = colors[rdg$winner]
n.sms = length(unique(rdg$Samples))


# make a region df that defines winner for each region
min_n_to_support = 5
rgn.df  = data.table(rdg %>% group_by(chr,start,end,name) %>%
                        summarise( CHM13 = sum( abs(CN-chm13) <= 1 ), GRCh38=sum(abs(CN-hg38) <= 1)) %>%
                        mutate( status = case_when( CHM13==0 & GRCh38==0 ~ "Neither", CHM13 >  GRCh38  ~ "CHM13", GRCh38 > CHM13 ~ "GRCh38", TRUE ~ "Polymorphic" )  ))
                        #mutate( status = case_when(CHM13 >= min_n_to_support & GRCh38 >= min_n_to_support ~ "Polymorphic", CHM13 >= min_n_to_support ~ "CHM13", GRCh38 >=min_n_to_support ~ "GRCh38", TRUE ~ "Neither" )  ))
rgn.df$region_color = colors[rgn.df$status]

# merge wiht super pop info and location
pop_location = data.table(super_pop = c("AMR","AFR","EA","SA","SIB","WEA","OCN"),
                          lat  = c( 13, 7.7, 33, 10.3, 61, 52, 1),
                          long  = c(-88, 21,  110, 75.5, 95, 27, 130)
); pop_location$lat = as.numeric(pop_location$lat); pop_location$long = as.numeric(pop_location$long)
rdg$sample = rdg$Samples
rdg=merge(rdg, popinfo, by="sample")
rdg= data.table( rdg %>% merge(pop_location, by="super_pop") %>% group_by(super_pop)%>% mutate(n_super_pop = length(unique(sample)))) 


# merge with region info
rdg = merge(rdg, rgn.df[,c("name","status")], on="name")
rdg$region_color = colors[rdg$status]

#
#
# This plot shows the cumlative distribution of CN imporvment
#
ecdf_rdg = pivot_longer(rdg, cols = which(colnames(rdg) %in% c("chm13_diff","hg38_diff")), names_to="ecdf_name", values_to = "ecdf_diff" )
ecdf_cn_plot = ggplot(data=ecdf_rdg) + 
  stat_ecdf(aes(x=abs(ecdf_diff), color=ecdf_name), size=3)+
  #[annotation_logticks()+
  #scale_x_continuous(trans="log10") +
  scale_color_manual(values=c(chm13_diff=NEWCOLOR, hg38_diff=OLDCOLOR))+
  coord_cartesian(xlim = c(0,30))+
  xlab("Maximum difference between sample CN and reference CN allowed to be considered correct")+
  ylab("Cumulative fraction of correct CN representation")+
  theme_cowplot()+theme(legend.position = "none")


#
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




pp=plotDefaultPlotParams(plot.type = 2); 
pp$plot.params$ideogramheight = 200
ideo = as.ggplot(expression(
  kp <- plotKaryotype(genome=GENOME, chromosomes = NOYM[NOYM %in% unique(rdg$chr)], cytobands = CYTO, plot.type = 2, plot.params = pp$plot.params),
  topp<-rdg[rdg$winner %in% c("CHM13") & abs(rdg$chm13_diff) <= 1 ,],
  botp<-rdg[rdg$winner %in% c("GRCh38") & abs(rdg$hg38_diff) <= 1,],
#kpPlotDensity(kp, data = toGRanges(rdg[winner=="Polymorphic",1:3]), col="gray")
kpPlotRegions(kp, data = toGRanges(topp[,c("chr","start","end")]), col=topp$color, ylim=n.sms, data.panel = 1, r0=.15),
kpPlotRegions(kp, data = toGRanges(botp[, c("chr","start","end")]), col=botp$color, ylim =n.sms, data.panel = 2, r0=.15),
kpPlotRegions(kp, data = toGRanges(rgn.df[rgn.df$status!="Polymorphic",c("chr","start","end")]), col=rgn.df$region_color[rgn.df$status!="Polymorphic"], ylim = n.sms, data.panel = "ideogram", r0=0.2, r1=0.8, border = NA)
)); ideo

TSpecial <- ttheme_minimal(
  colhead=list(fg_params=list(fontface=4L)),
  rowhead=list(fg_params=list(fontface=3L)))

t1 = rdg %>% group_by(winner) %>% summarise(count=length(winner))
t2 = rgn.df %>% group_by(status) %>% summarise(count = length(status))

p3 = ggplot(data=t1) + geom_bar(aes(y=winner, x=count, fill=winner), stat = "identity") + 
  geom_text(aes(label=count, y=winner, x=count), vjust=-1, angle=-90) +
  scale_fill_manual(values=colors) + theme_cowplot() + xlab("# of sample genotypes with CN closer to CHM13 or GRCh38") + ylab("") + theme(legend.position = "none") ;p3
p4 = ggplot(data=t2) + geom_bar(aes(y=status, x=count, fill=status), stat = "identity") + 
  geom_text(aes(label=count, y=status, x=count), vjust=-1, angle=-90) +
  scale_fill_manual(values=colors) + theme_cowplot() + xlab(glue("# of regions where more samples genotype exactly with either CHM13 or GRCh38")) + ylab("") + theme(legend.position = "none"); p4




#
# world map
#
sizescale= 2
zz_extra = data.table(rdg %>% group_by(super_pop,lat,long,n_super_pop,chr,start,end,name) %>%
                        summarise( CHM13 = sum( abs(CN-chm13) <= 1 ), GRCh38=sum(abs(CN-hg38) <= 1)) %>%
                        mutate( status = case_when( CHM13==0 & GRCh38==0 ~ "Neither", CHM13 >  GRCh38  ~ "CHM13", GRCh38 > CHM13 ~ "GRCh38", TRUE ~ "Polymorphic" )) %>%
                        ungroup() %>% group_by(super_pop,lat,long,n_super_pop) %>%
                        summarise(CHM13 = sum(status=="CHM13")/length(status), GRCh38 = sum(status=="GRCh38")/length(status),Polymorphic = sum(status=="Polymorphic")/length(status), radius=sizescale*sqrt(unique(n_super_pop))));#unique(n_super_pop)/sizescale));zz_extra



world <- map_data('world')
map <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill=NA, color="black") +
  coord_quickmap() +
  geom_scatterpie(aes(x=long, y=lat, group=super_pop, r=radius),
                     data=zz_extra, cols=c("CHM13", "GRCh38", "Polymorphic"), color=NA, alpha=0.95) +
  geom_scatterpie_legend(zz_extra$radius, x=-160, y=-55, labeller = function(x){x*x/sizescale} ) +
  scale_fill_manual(values=colors) + theme_map() + ylab("Latitude") + xlab("longitude") + theme(legend.position = "top", legend.title = element_blank(), legend.justification = "center"); map


#
# mbp histogram
# 
bp_hist = ggplot(data=rgn.df,aes(x=status, weight=end-start, fill=status))+
  geom_bar()+
  theme_cowplot()+
  scale_y_continuous(labels = comma) + 
  scale_fill_manual(values=colors) +
  geom_text(stat='count', aes(label=comma(..count..)), vjust=-1)+
  xlab("") + ylab("# of bp where more samples genotype exactly") +
  theme(legend.position = "none")

#
# Final figure
#
fig_top = plot_grid(
  plot_grid(p4, ideo, rel_heights = c(.4,1), nrow=2),
  plot_grid(p3, p,    rel_heights = c(.4,1), nrow=2),
  rel_widths = c(1,1), labels = c("a","b"))
fig_bot = plot_grid(bp_hist, ecdf_cn_plot, rel_widths = c(1,2), labels =c("c","d"))
fig = plot_grid(fig_top, fig_bot, ncol=1,labels = c(NA,NA), rel_heights = c(4,3))

# simplified figure
fig = plot_grid(plot_grid(bp_hist,p, labels=c("a","b")), ecdf_cn_plot, ncol = 1, labels =c(NA,"c"))


scale=0.7
ggsave(glue("figures/wssd_{pname}.pdf"), plot=fig, height = 16*scale, width = 16*scale)
ggsave(glue("figures/wssd_{pname}.png"), plot=fig, height = 16*scale, width = 16*scale, dpi=300)
fig
