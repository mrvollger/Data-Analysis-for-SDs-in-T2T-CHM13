#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")



#bed = fread('../sd_regions_in_hifi_wga/lpa/all.gene.models.bed', col.names = c("contig", "starto","endo", "target", "qual", "strand","sample","hap"));bed$len = bed$endo - bed$starto; bed

cnbed=fread("lpa_figure/LPA.aln.tbl")
cnbed$contig = cnbed[["reference_name"]]

cnbed = data.table(cnbed %>% group_by(contig) %>% summarise(CN = n()))
#cnbed$sample = gsub("\\..*","", cnbed$contig)
#cnbed$hap = gsub("__.*","",gsub(".*\\.", "", cnbed$contig))
#cnbed

bed =fread("lpa_figure/LPA.bed", col.names = c("contig", "starto", "endo"))
bed = merge(bed, cnbed, all.x=TRUE, by="contig")


bed$sample = gsub("__.*", "" , gsub("LPA__", "", bed$contig))
bed$hap=gsub("__.+", "", gsub("LPA__[^_]+__", "", bed$contig))





bed = merge(bed, fread("../assemblies_for_anlysis/sample_info/Master.tbl"), all.x=T, by=c("sample","hap")); bed
bed[hap=="pri"]$SuperPop=""
bed[hap=="pri"]$hap=""
bed$type="Human"

bed[grepl("_",bed$sample)]$type = "NHP"
bed[grepl("_",bed$sample) & grepl("alt", bed$contig) ]$hap = "alt"
bed[grepl("_",bed$sample) & grepl("pri", bed$contig) ]$hap = "pri"
#  
#gsub(".*__", "",tochange)

bed$lab = gsub("_"," ",paste(bed$sample, bed$hap))#, bed$Population, sep="_")


bed = data.table(bed %>% group_by(lab) %>% mutate( clength = max(endo)))
bed$starto2 = bed$clength - bed$endo
bed$endo2 = bed$clength - bed$starto

bed[bed$sample == "HG002",]$SuperPop = "EUR"
bed = data.table(bed %>% group_by(lab) %>% mutate( clength = max(endo2), start = starto2 - min(starto2), end = endo2-min(starto2)))

at_risk=data.table(bed %>% group_by(sample,hap) %>%
  summarise(at_risk = max(CN) <=20 ) )


bed=merge(bed, at_risk, by=c("sample","hap"), all.x=T)
bed$at_risk = as.character(bed$at_risk)
bed[type=="NHP"]$at_risk = "Unknown (NHP)"


bed$SuperPop = factor(bed$SuperPop, levels = sort(unique(bed$SuperPop)))

bed$sort = bed$clength
bed[sample == "CHM13"]$sort = 0
bed[sample == "GRCh38chrOnly"]$sort = 1
#bed[sample == "HG03125"]$sort = 1e6

bed = bed[order(type,sort), ]
bed$lab = factor(bed$lab, levels = rev(unique(bed$lab)))



#bed %>% group_by(sample) %>% mutate(Dip_CN = sum(unique(CN)), nhaps = length(unique(hap)))

p=ggplot(data=bed, aes(color=at_risk) )+
  geom_segment(aes(x=start, xend=end, y=lab, yend=lab), size=6)+
  #scale_x_continuous( expand = c(0, 0), limits = c(1,max(df$end)),breaks = c(1,max(df$end)) ) +
  geom_segment(data = bed %>% group_by(lab,SuperPop) %>% summarise(start=min(start), end=max(end), at_risk=unique(at_risk)), aes(x=start, xend=end, y=lab, yend=lab), size=0.75)+
  ylab('')+xlab('')+
  theme_classic()+
  geom_text(data = unique(bed[,c("lab","CN","at_risk")]), aes(x=-1e4, y=lab, label=CN))+
  geom_text(aes(x=-1e4, y=max(as.numeric(lab))+1, label="CN"), color="black")+
  scale_x_continuous(label = comma)   +
  scale_color_manual(values = c(`TRUE`=NEWCOLOR, `FALSE`="black", `Unknown (NHP)`="#4B0082"))+
  #scale_color_brewer(palette = "Dark2")+
  #facet_grid(lab~.)+
  #scale_y_discrete(limits = c(0, max(length(unique(bed$lab)))) )+
  theme_cowplot()+
  coord_cartesian(clip = 'off')+
  theme(plot.margin=margin(1,1,1,1, 'cm'), axis.line.y =  element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")+
  guides(color=guide_legend(title ="At increased risk for CVD"));p
ggsave( "lpa_figure/lpa_gene_models.pdf", plot=p, height = 9, width = 16)



#g=plot_grid(p2, NULL, p, labels = c('a','c', 'b'));g
#ggsave("fig.pdf", plot = g, height = 12, width = 16)

#ggdraw() + draw_image("aln.png", scale = 0.7)



#paf = fread("../minigraph/lpa/out.paf", fill=T, header=F)[,1:12]; colnames(paf) = c("q", "qle", "qst", "qen", "strand", "t", "tle", "tst", "ten", "x","y","z"); paf; unique(paf$t)
#ggplot(data=paf) +geom_segment(aes(x=tst, xend=ten, y=q, yend=q, color=q)) + facet_wrap(vars(t), ncol=length(unique(paf$t)), scales="free_x") +theme_cowplot() + theme(legend.position = "none")



