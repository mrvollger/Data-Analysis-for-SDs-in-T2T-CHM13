#!/usr/bin/env Rscript

duptaged = rgntag(DM, SEDEF, "SD")

#
#
#
#

dups = data.table(rbind(DM, DM_38) %>% 
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
p2 = ggplot(data = dups[1:n], aes(y=label, weight=count/1000, fill=Assembly)) + 
  geom_bar(position="dodge") +
  scale_fill_manual(values=COLORS)+
  xlab("Total duplicon content (kbp)") +
  ylab("Duplicon (Gene:count,...)") +
  scale_y_discrete(position = "right")+
  scale_x_reverse(label=comma)+
  labs(title="30 duplicons with largest change between\nT2T CHM13 and GRCh38") +
  theme_cowplot()+theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5));p2



#
#
#
#

t =bind_rows(DupMasker=duptaged[,1:3],.id="Type",
  `SDs without DupMasker`=sedef[!sedef$Duplicon, 1:3],
  `30 most expanded duplicons`=duptaged[duptaged$duplicon %in% dups[1:n]$duplicon, 1:3] )
pal=c(`DupMasker`="darkgrey",`SDs without DupMasker`="red",`30 most expanded duplicons`="darkgreen")
p=ggplot()+
  #geom_density(data=t[chr!="chrMT"], aes(x=start/1000000, weight=(end-start), fill=Type, group=Type), alpha=0.5,adjust = 0.5)+
  geom_histogram(data=t[chr!="chrMT"], aes(x=start/1000000, weight=(end-start), fill=Type, group=Type), 
                 binwidth = 2, alpha=0.75, position="stack")+
  #bgeom_segment(data=FAI, aes(x=0,xend=max(chrlen)/1000000,y=0,yend=0,group=chr))+
  scale_x_continuous(labels = comma)+
  scale_fill_manual(values=pal)+
  facet_wrap(chr~., ncol=2)+
  ggtitle("SDs unannotated by DupMakser (Red)") +
  xlab("Genomic position (Mbp)")+
  theme_cowplot()+theme(legend.position = "bottom")+
  theme(strip.text = element_text(size = 10, margin = margin()), 
        strip.background = element_blank(),
        legend.title = element_blank())
p


#
#
#
#

pal = dupmasker_pal()[]
pal = c(pal, "ALR/Alpha"="#e50000", "BSR/Beta"="#4c4cff", "SAR"="#555555", "TRUE"="#000000", "FALSE"="#555555", 
        "HSATII"="#4c4cff", "CER"="#555555", "GSATII"="#4c4cff") 
achro_rm = achro(RM)[type=="Satellite" & (end-start)>0 & t1 %in% names(pal)] 

labels= data.table(x=0, y=c(.5, .15,.05, -.15), label=c("DupMasker","SDs","Ancient SDs","Satellites"))
p3 = ggplot()+
  geom_rect(data=achro_rm, aes(xmin=start, xmax=end, ymin=-.25, ymax=-0.05, fill=t1))+
  geom_rect(data=achro(DM), aes(xmin=start, xmax=end, ymin=0.25, ymax=1, fill=duplicon))+
  geom_rect(data=achro(SEDEF), aes(xmin=start, xmax=end, ymin=0.1, ymax=0.2, fill=intra))+
  geom_rect(data=achro(LOW), aes(xmin=start, xmax=end, ymin=0, ymax=0.1), fill="lightgreen", alpha=0.5)+
  geom_text(data=labels, aes(x=x,y=y,label=label) , hjust="right")+
  facet_wrap(chr~., ncol=1)+
  scale_fill_manual(values=pal, labels = names(pal))+
  theme_cowplot()+theme(legend.position = "none")+
  xlab("Genomic position")+scale_x_continuous(labels = comma)+
  ylab("")+
  theme(strip.text = element_text(size = 16, margin = margin()), strip.background = element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y = element_blank()
        ); p3
  
top = plot_grid(p2,p, ncol=2, labels = c("a","b"))+theme(plot.title = element_text(hjust = 0.5))
x = plot_grid(top, p3, nrow=2,labels = c("","c"))
ggsave("duplicons2.pdf", plot=x, width=24, height = 16)

x = plot_grid(p2,p,p3, nrow=3,labels = c("a","b","c"))
ggsave("duplicons_long.pdf", plot=x, width=24, height = 20)


