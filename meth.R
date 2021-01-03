#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")
df = fread("../t2t_globus_share/team-epigenetics/20200727_methylation/v1.0_methylation/t2t_chm13v1.0_SD_clustered_methylation.bed")

mycol = function(x){
  if(x["clust_meth"]=="meth"){
    val = 1-as.numeric( x["avgmeth"] )
    return(transparent(NEWCOLOR,val))
  }else{
    val = as.numeric( x["avgmeth"] )
    return(transparent(OLDCOLOR, val))
  }
}

df$color = apply(df,1, mycol)
meth<-toGRanges(df[clust_meth=="meth"])
unmeth<-toGRanges(df[clust_meth!="meth"])


b=as.ggplot(expression(
  kp <- plotKaryotype(genome = GENOME, chromosomes =NOYM, cytobands = CYTO, plot.type = 6),
  kpPlotRegions(kp, data=meth, col=meth$color, border = NEWCOLOR),
  kpPlotRegions(kp, data=unmeth, col=unmeth$color, border = OLDCOLOR)
))+theme_nothing()

a= ggplot(data=df) + 
  geom_histogram(aes(avgmeth, fill=color, weight=end-start), bins=50) + 
  xlab("Average methylation")+ylab("# of methylated bp") +
  scale_y_continuous(labels = comma)+
  scale_fill_identity()+
  theme_cowplot()

l=get_legend(ggplot(data=data.frame(Status=c("Centromere","Methylated","Unmethylated")))+geom_bar(aes(0,fill=Status)) + 
               scale_fill_manual(values = c("Black",NEWCOLOR,OLDCOLOR))+ theme_cowplot() + theme(legend.position = "top", legend.justification = "center", legend.title = element_blank()))

meth_block_plot = plot_grid(a,b,l,ncol = 1, rel_heights = c(1.5,3,.25), rel_widths = c(.8,1,1)) 
scale=1.25
ggsave("supp/meth_locations.pdf", plot=meth_block_plot, height = 8*scale, width = 12*scale)





#
#
# meth genes plots
#
#
in.df = METH_SD_GENES[, c("width","num_motifs_in_group","called_sites", "group_sequence", "total_calls","called_sites_methylated"):=NULL]
# summarize data by gene
in.df$quartile = cut(in.df$n_transcripts, breaks = c(0,1,10,Inf), include.lowest = T, right = F)
n_windows=250
in.df$cut = cut(in.df$dist, breaks=n_windows)
in.df$is_sd = in.df$is_sd > 0
grouped.df = in.df %>%
  group_by(gene, n_transcripts, quartile, cut, is_sd, genewidth) 

gene.df = grouped.df %>%
  summarise(methylated_frequency = median(methylated_frequency), 
            n_cpg = n(), 
            chr=chr[1], 
            start=min(start), 
            end=max(end)) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% # pull out the start and end of the cut intervals
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
gene.df

# summarize data along all genes
gene.ave.df = gene.df %>% ungroup() %>% 
  group_by(is_sd, cut, quartile,min,max) %>%
  summarise(
    med = median(methylated_frequency), 
    top = quantile(methylated_frequency, 0.75),
    bot = quantile(methylated_frequency, 0.25),
    med_cpg = mean(n_cpg/(end-start+1))) 

# count num of genes
n_gene = gene.df %>% 
  group_by(quartile, is_sd) %>%
  summarise(n_genes=length(unique(gene)))

# main figure
iso_meth_plot = ggplot(data=gene.ave.df, 
       aes(x=min, color=is_sd, fill=is_sd)) +
  # add the data intervals
  geom_ribbon(aes(ymin=bot, ymax=top), color = transparent("black",1), alpha=0.2, size=0)+
  # plot the methylation freq
  geom_point(aes(y=med), fill="black", alpha=1, size=1) +
  geom_line( aes(y=med), fill="black",alpha=1, size=.5)+
  # add the cpg dentity
  #geom_line(aes(y=med_cpg/max(gene.ave.df$med_cpg)), alpha=1, size=.1)+
  geom_text(data = n_gene,
            aes(x=c(2-is_sd),
                y=c(1),
                label=paste("# genes =",comma(n_genes))
                ), 
            hjust=0)+
  geom_vline(xintercept=c(0,2), linetype=5 )+
  scale_color_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  scale_fill_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  scale_x_continuous(breaks = c("-10 kbp"=-1, "TSS"=0, "TTS"=2, "+10 kbp"=3)   )+
  facet_wrap(vars(quartile), nrow = 3)+
  theme_cowplot() +
  theme(legend.position = "bottom") +
  ylab("CpG methylation frequency")+
  xlab("Normalized position along gene body")+
  labs(subtitle="# Iso-Seq transcripts [min, max)"); iso_meth_plot
  
# controlling for cpg density
iso_ncpg_plot = ggplot(data=gene.ave.df, 
                       aes(x=min, color=is_sd, fill=is_sd)) +
  geom_text(data = n_gene,
            aes(x=c(2-is_sd),
                y=c(1.1),
                label=paste("# genes =",comma(n_genes))
            ), 
            hjust=0)+
  # plot the methylation freq
  geom_point(aes(y=med_cpg), alpha=1, size=1) +
  geom_line( aes(y=med_cpg), alpha=1, size=.5)+
  # add the cpg dentity
  #geom_line(aes(y=med_cpg/max(gene.ave.df$med_cpg)), alpha=1, size=.1)+
  geom_vline(xintercept=c(0,2), linetype=5 )+
  scale_color_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  scale_fill_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  scale_x_continuous(breaks = c("-10 kbp"=-1, "TSS"=0, "TTS"=2, "+10 kbp"=3)  )+
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  facet_wrap(vars(quartile), ncol = 3)+
  theme_minimal_hgrid() +
  theme(legend.position = "none") +
  ylab("Average density of CpG sites")+
  xlab("Normalized position along gene body")+
  labs(subtitle="# Iso-Seq transcripts [min, max)"); iso_ncpg_plot

# calculate the legnths of the genes
tmp.df = unique(gene.df[,c("genewidth","gene", "is_sd", "quartile")])
length_plot = ggplot(data = tmp.df,
                     aes(x=genewidth,y=is_sd, fill = is_sd, color=is_sd)
                     ) +
  #geom_density_ridges(color="black")+
  #geom_rug()+
  geom_histogram(aes(genewidth, y=..density.., fill=is_sd), bins=30, alpha=0.7, position = "identity")+
  scale_color_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  scale_fill_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  scale_x_continuous(label=comma, trans="log10")+
  annotation_logticks(sides = "b")+
  ylab("Density")+ xlab("Gene length") +
  facet_grid(.~quartile)+
  #facet_zoom(x = genewidth > 10000 & genewidth < 100000)+
  theme_minimal_vgrid()+
  theme(legend.position = "bottom"); length_plot

ggsave("supp/meth_ave_num_cp.pdf", height = 12, width = 20, plot=plot_grid(iso_ncpg_plot, length_plot, nrow=2))

sd_genes = unique(in.df[in.df$is_sd], by= c("chr","gene","quartile", "n_transcripts"))
sd_genes$origin = gsub("_.*","",sd_genes$gene)
sd_genes_0 = sd_genes[n_transcripts==0]
sd_genes_0$prefix = substr(sd_genes_0$origin,1,4) #gsub( '(.*[[:alpha:]])([0-9]{1,3})$', '\\1',  sd_genes_0$origin)
sd_genes_0 = data.table(sd_genes_0 %>% group_by(prefix) %>% arrange(origin) %>% summarise(count = n(), gene=origin[1], start=start[1], chr=chr[1]))

sd_genes_0$y = runif(nrow(sd_genes_0),.3,.7)
sd_genes_0$p = 1
sd_genes_0$p[seq_along(sd_genes_0$p) %% 2 > 0] = 2
d1 = sd_genes_0[p==1]
d2 = sd_genes_0[p==2]

make_gene_plot <- function(chrs){
  tchrs <<- chrs
  m= c(0.1,.1,0.5)
  genes = as.ggplot(expression(
    kp <- plotKaryotype(genome=GENOME, cytobands = CYTO, chromosomes = tchrs, plot.type = 3),
    kpPlotMarkers(kp, chr = d1$chr, x = d1$start, labels = d1$gene, y=d1$y, adjust.label.position = T, marker.parts = m, data.panel = 1, ignore.chromosome.ends = T, label.dist=0, label.margin = 0, cex=0.75),
    kpPlotMarkers(kp, chr = d2$chr, x = d2$start, labels = d2$gene, y=d2$y, adjust.label.position = T, marker.parts = m, data.panel = 2, ignore.chromosome.ends = T, label.dist=0, label.margin = 0, cex=0.75)
  ))
  genes
}
#genes=make_gene_plot(NOYM[1:9])
#genes2=make_gene_plot(NOYM[10:23])
#plot_grid(iso_meth_plot, genes,genes2, rel_heights = c(2,2,2), ncol=1)


gage = gene.df[grepl("TBC1D3([[:alpha:]])*$", gene.df$gene),] 
gage = gene.df[grepl("NPIPA", gene.df$gene),] 
gene_plot = ggplot(data=gage %>% arrange(n_transcripts))+
  geom_vline(xintercept=c(0,2), linetype=5 )+
  geom_point(aes(x=min, y=methylated_frequency), size=0.5, alpha=0.5)+
  geom_line(aes(x=min, y=rollmean(methylated_frequency, 5, na.pad=TRUE)), size=.5, color=NEWCOLOR)+
  facet_wrap(n_transcripts~gene, nrow=1)+
  scale_x_continuous(breaks = c("TSS"=0, "TTS"=2)   )+
  theme_cowplot() +
  theme(legend.position = "none") +
  ylab("CpG methylated frequency")+
  xlab("Normalized position along gene body");gene_plot

#meth_gene_plot = plot_grid(iso_meth_plot, gene_plot, rel_heights = c(2,2), ncol=1)
meth_fig =plot_grid(plot_grid(meth_block_plot, iso_meth_plot, rel_widths = c(1,2), labels = c("a","b")), gene_plot, nrow=2, rel_heights = c(2,1), labels = c(NA,"c"))
meth_fig
ggsave("figures/meth_fig.pdf", plot=meth_fig, height = 12, width = 16)
