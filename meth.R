#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")


df=METH_CLUSTERS




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
  kp <- plotKaryotype(genome = GENOME, chromosomes =NOYM, cytobands = CYTO[CYTO$gieStain != "stalk"], plot.type = 6),
  kpPlotRegions(kp, data=meth, col=meth$color, border = NEWCOLOR),
  kpPlotRegions(kp, data=unmeth, col=unmeth$color, border = OLDCOLOR)
))+theme_nothing(); b

a= ggplot(data=df) + 
  geom_histogram(aes(avgmeth, fill=color, weight=end-start), bins=50) + 
  xlab("Average methylation")+ylab("# of methylated bp") +
  scale_y_continuous(labels = comma)+
  scale_fill_identity()+
  theme_cowplot()

leg = c("Methylated","Unmethylated","Centromere")#,"Gap")#,"#647FA4"
l=get_legend(
  ggplot(
    data=data.frame( Status=factor(leg, levels = leg) ) 
    ) + 
    geom_bar(aes(0,fill=Status)) + 
    scale_fill_manual(values = c(NEWCOLOR,OLDCOLOR,"black"))+ theme_cowplot() + theme(legend.position = "top", legend.justification = "center", legend.title = element_blank()))

meth_block_plot = cowplot::plot_grid(l,b,a, ncol = 1, rel_heights = c(0.25,3,1.5)) 
meth_block_plot
scale=1.25

ggsave(glue("{SUPP}/meth_ideogram.pdf"), plot=meth_block_plot, height = 8*scale, width = 12*scale)


#
#
#
# meth genes plots
#
#
#
SAMPLE="HG002"
SAMPLE=""
if(SAMPLE=="HG002"){
  meth_and_transcript = METH_SD_GENES_002
}else{
  meth_and_transcript = METH_SD_GENES
}
in.df = meth_and_transcript

# summarize data by gene
in.df$quartile = cut(in.df$n_transcripts, breaks = c(0,1,10,Inf), include.lowest = T, right = F)
n_windows=200
in.df$cut = cut(in.df$dist, breaks=n_windows)
in.df$is_sd = in.df$is_sd > 0
cluster <- new_cluster(12)
grouped.df = in.df %>%
  group_by(gene, n_transcripts, quartile, cut, is_sd, genewidth) 
#%>% 
 # partition(cluster)

gene.df = grouped.df %>%
  summarise(methylated_frequency = median(methylated_frequency), 
            n_cpg = n(), # might need to change back to n()
            chr=chr[1], 
            start=min(start), 
            end=max(end)) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% # pull out the start and end of the cut intervals
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) 
#%>%
 # collect()


gene.df
unique(gene.df[,c("chr","start","end","cut","gene")])
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
       aes(x=max, color=is_sd, fill=is_sd)) +
  # plot the methylation freq
  geom_point(aes(y=med), alpha=1, size=1,fill="black") +
  geom_line( aes(y=med), alpha=1, size=.5)+
  # add the cpg dentity
  #geom_line(aes(y=med_cpg/max(gene.ave.df$med_cpg)), alpha=1, size=.1)+
  geom_text(data = n_gene,
            aes(x=c(2-is_sd+.1),
                y=c(0.05),
                label=paste("# genes =",comma(n_genes))
                ),
            size=5,
            show.legend = FALSE,
            hjust=0)+
  # add the data intervals
  geom_ribbon(aes(ymin=bot, ymax=top), color = transparent("black",1), alpha=0.2, size=0)+
  
  geom_vline(xintercept=c(0,2), linetype=5)+
  scale_color_manual("Median\nmethylation", values = c(OLDCOLOR, NEWCOLOR), labels=c("unique", "SD"))+
  scale_fill_manual("Methylation\nquartiles (1-3)", values = c(OLDCOLOR, NEWCOLOR), labels=c("unique", "SD"))+
  scale_x_continuous(breaks = c("-10 kbp"=-1, "TSS"=0, "TTS"=2, "+10 kbp"=3)   )+
  
  guides(color=guide_legend(order=1))+
  
  facet_wrap(vars(quartile), nrow = 3)+
  theme_cowplot() +
  theme(legend.position = "right", legend.justification = "center") +
  ylab("CpG methylation frequency")+
  xlab("Normalized position along gene body")+
  labs(subtitle="# Iso-Seq transcripts [min, max)"); iso_meth_plot
  
#
#
#
# controlling for cpg density
#
#
#
#just_genes = readbed("data/methylation_sd_analysis/sd.transcripts.and.meth.bed","z")
just_genes = METH_JUST_GENES

strand_scale = function(just_genes,  width = 0.01, scale = 1, offset=0){
  mmin = 0
  mmax = 1
  data.table(tidyr::crossing(just_genes, dist_min = seq(mmin, mmax - width, width) ) %>%
    mutate(dist_min = case_when( strand == "+" ~ dist_min,
                                 strand == "-" ~  mmax - dist_min - width),
           dist_max = dist_min + width,
           length = end - start
           ) %>% 
    mutate(start2 = format(ceiling(start + length * dist_min), scientific = FALSE, trim=TRUE),
           end2 = format(ceiling(start + length * dist_max), scientific=FALSE, trim=TRUE) 
           ) %>% 
    mutate(start = start2, end = end2,
           dist_min = dist_min * scale + offset,
           dist_max = dist_max * scale + offset) %>%
    filter(start >= 0) %>%
    dplyr::select(-c(start2, end2)))
}
left_slop <- function(bed, slop=10000) {
  df = copy(bed)
  df$start = bed$start - slop
  df$end = bed$start
  df
}
right_slop <- function(bed, slop=10000) {
  df = copy(bed)
  df$end = bed$end + slop
  df$start = bed$end
  df
}
up_slop <- function(bed){
  rbind(left_slop(bed[strand=="+"]), right_slop(bed[strand=="-"]))
}
down_slop <- function(bed){
  rbind(right_slop(bed[strand=="+"]), left_slop(bed[strand=="-"]))
}

scaled_genes = rbind(strand_scale(up_slop(just_genes), width = 0.02, offset=-1),
      strand_scale(down_slop(just_genes), width = 0.02, offset=2),
      strand_scale(just_genes, scale=2))[order(chr,start)]

if(F){
  write.table(scaled_genes, file="data/scaled_binned_slop_genes.bed", 
              col.names = F, sep="\t",row.names = F, quote = F)
  FA="../Assembly_analysis/SEDEF/chm13.draft_v1.0_plus38Y_masked.fasta"
  system(glue("bedtools nuc -fi {FA} -C -pattern CG -bed data/scaled_binned_slop_genes.bed > data/scaled_binned_slop_genes.cpg.bed"))
}

cpg_data = CPG_DATA
colnames(cpg_data)[1:length(colnames(scaled_genes))] = colnames(scaled_genes) 
cpg_data$n_cpg = cpg_data$`27_user_patt_count`
cpg_data$quartile = cut(cpg_data$n_transcripts, breaks = c(0,1,10,Inf), include.lowest = T, right = F)
cpg_data$is_sd = cpg_data$is_sd > 0


cpg_density = cpg_data %>% 
  group_by(is_sd, quartile, dist_min, dist_max) %>%
  #filter(  ) %>%
  summarise(med_cpg = mean(n_cpg/(end-start)), n_cpg = median(n_cpg)) 

iso_ncpg_plot = ggplot(data=cpg_density, 
                       aes(x=dist_max, color=is_sd, fill=is_sd)) +
  geom_text(data = n_gene,
            aes(x=c(2-is_sd),
                y=c(0),
                label=paste("# genes =",comma(n_genes))
            ), 
            hjust=0)+
  # plot the methylation freq
  geom_point(aes(y=med_cpg), alpha=1, size=1) +
  geom_line( aes(y=med_cpg), alpha=1, size=.5)+
  #geom_smooth(aes(y=med_cpg))+
  # add the cpg dentity
  #geom_line(aes(y=med_cpg/max(gene.ave.df$med_cpg)), alpha=1, size=.1)+
  geom_vline(xintercept=c(0,2), linetype=5 )+
  scale_color_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  scale_fill_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  scale_x_continuous(breaks = c("-10 kbp"=-1, "TSS"=0, "TTS"=2, "+10 kbp"=3)  )+
  #scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  facet_wrap(vars(quartile), ncol = 1)+
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

ggsave(glue("{SUPP}/meth_ave_num_cp{SAMPLE}.pdf"), height = 12, width = 20, plot=cowplot::plot_grid(iso_ncpg_plot, length_plot, nrow=2, rel_heights = c(2,1)))

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

n_rolling_mean = 10
gage = gene.df[grepl("TBC1D3([[:alpha:]])*$", gene.df$gene),] 
gage = gene.df[grepl("NPIPA", gene.df$gene),] 
#gage = gene.df[grepl("NBPF", gene.df$gene),] 

gene_plot = ggplot(data=gage %>% arrange(n_transcripts))+
  geom_vline(xintercept=c(0,2), linetype=5 )+
  geom_point(aes(x=min, y=methylated_frequency, fill=""), size=0.5, alpha=1)+
  scale_fill_manual("Data", values = "black")+
  geom_line(aes(x=min, y=rollmean(methylated_frequency, n_rolling_mean, na.pad=TRUE, align = "center"), color=""), size=.5)+
  scale_color_manual(glue("Rolling mean (n = {n_rolling_mean})"), values = NEWCOLOR)+
  facet_wrap(n_transcripts~gene, nrow=1)+
  scale_x_continuous(breaks = c("TSS"=0, "TTS"=2)   )+
  theme_cowplot() +
  theme(legend.position = "top") +
  ylab("CpG methylation frequency")+
  xlab("Normalized position along gene body");gene_plot

#meth_gene_plot = plot_grid(iso_meth_plot, gene_plot, rel_heights = c(2,2), ncol=1)
meth_fig = cowplot::plot_grid(
  cowplot::plot_grid(
    meth_block_plot, iso_meth_plot,
    rel_widths = c(1,2), labels = c("a","b")
    ), 
  gene_plot, nrow=2, rel_heights = c(2,1), labels = c(NA,"c"))

meth_fig
ggsave(glue("{SUPP}/meth_fig{SAMPLE}.pdf"), plot=meth_fig, height = 12, width = 16)




#
# near tts site methylation
#
near_tss  = gene.df %>% 
  group_by(is_sd, quartile) %>%
  mutate(median = median(methylated_frequency)) %>%
  filter(max >= -0.04 & max <= 0.04); near_tss

sdtts = near_tss[near_tss$is_sd & near_tss$n_transcripts == 0,]
utts = near_tss[!near_tss$is_sd & near_tss$n_transcripts == 0,]

ggplot(data=near_tss,
       aes(methylated_frequency,
           color=is_sd, fill=is_sd)) +
  geom_density() +
  geom_histogram(bins=30)+
  facet_grid(is_sd~quartile) +
  #geom_line(aes(y=cut))+
  scale_color_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  scale_fill_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  theme_cowplot()
  


test_hyper_methylation = wilcox.test(sdtts$methylated_frequency, 
            utts$methylated_frequency, alternative = "greater")

cat(as.character(test_hyper_methylation), file=glue("{SUPP}/meth_hyper_over_tss{SAMPLE}.txt"), sep="\n")

z = near_tss %>% group_by(is_sd, quartile) %>% 
  summarise(median_at_tts = median(methylated_frequency), median_overall = unique(median)) %>%
  mutate(difference = median_at_tts - median_overall)
z
write.table(z, file=glue("{SUPP}/meth_hyper_over_tss{SAMPLE}.tbl"))






#
#
#
# look for cpg islands
#
#
#
gene_tbl = unique(data.table(gene.df[,c("gene","n_transcripts","quartile","is_sd")]))
bins=1500
x=fread(glue("data/flanks.{bins}.gc.content.bed"))
genes_with_gc = merge(gene_tbl, x, by.x=c("gene"), by.y=c("4_usercol") )

genes_with_gc$`CpG count` = genes_with_gc$`17_user_patt_count`
genes_with_gc$`GC content` = genes_with_gc$`9_pct_gc`
genes_with_gc = data.table(genes_with_gc %>% pivot_longer(c("CpG count", "GC content")))
t = copy(genes_with_gc)
t$quartile = "All"
genes_with_gc = rbind(genes_with_gc, t)


stats =   genes_with_gc %>% group_by(is_sd , quartile, name) %>%
  filter( name == "GC content") %>%
  summarise(mean=mean(value), median = median(value), name="GC content")
stats

gc_plot = ggplot(data=genes_with_gc, aes(x=value, fill=is_sd, color=is_sd))+
  #geom_histogram(aes(y=..density..), position = "identity", bins=min(bins+1,51), alpha=0.5)+
  geom_density_ridges(aes(y=is_sd), alpha=0.8)+
  geom_count(aes(y=is_sd))+
  #scale_size()+
  geom_vline(data=stats, aes(xintercept=median,color=is_sd), size=2, alpha=0.8)+
  geom_vline(data=data.frame(x=0.41, name="GC content"), aes(xintercept=x), size=2, color="gray")+
  facet_grid(quartile~name, scales = "free")+
  theme_cowplot()+
  theme(legend.position = "bottom")+
  scale_color_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  scale_fill_manual("SD gene", values = c(OLDCOLOR, NEWCOLOR))+
  xlab(glue("CpG count (left) and GC content (right) {bins} bp upstream and downstream of TSS"))+
  ylab("Binned Iso-Seq transcript counts")
gc_plot 
ggsave(glue("{SUPP}/meth_tss_gc_content{SAMPLE}.pdf"), plot=gc_plot, height = 12, width = 16)





#
#
#
# check if the untranscribed genes overlap with hypomethylated blocks
#
#
#
gene_meth_pos = unique(in.df[,c("chr","gene_start","gene_end","gene", "n_transcripts", "is_sd")])
test_un = gene_meth_pos[is_sd & n_transcripts==0 & gene %in% near_tss$gene]
test_un

unmeth_x = unmeth[unmeth$avgmeth < 0.30]

size_meth = sum(width(GenomicRanges::reduce(meth))); size_meth
size_unmeth = sum(width(GenomicRanges::reduce(unmeth_x))); size_unmeth
size_total = size_meth + size_unmeth

inter_unmeth = subsetByOverlaps(toGRanges(test_un), toGRanges(unmeth_x))
inter_meth = subsetByOverlaps(toGRanges(test_un), toGRanges(meth))


rand_n_not_meth = function(fake){
  results = sample(c(0, 1), size=length(inter_meth)+length(inter_unmeth), replace=T, prob=c(size_meth/size_total, size_unmeth/size_total))
  sum(results)
}
permutation_results = sapply(seq(10000), rand_n_not_meth)
our_result = length(inter_unmeth)
hist(permutation_results, breaks=50)#, breaks=seq(125,225))
lines(x = rep(our_result, 2), y=c(0,1000000), col="red")

p_value_for_unmth_enrichment = 1 - sum(our_result > permutation_results)/ length(permutation_results)
p_value_for_unmth_enrichment

cat(p_value_for_unmth_enrichment,file=glue("{SUPP}/meth_hypo_block_pvalue{SAMPLE}.txt"),sep="\n")





































if(FALSE){
#
# motif ?
#
motif_df = meth_and_transcript[group_sequence != "split-group" & dist > -.002 & dist < .002]
motif_df$A = str_count(motif_df$group_sequence,"A")
motif_df$T = str_count(motif_df$group_sequence,"T")
motif_df$C = str_count(motif_df$group_sequence,"C")
motif_df$G = str_count(motif_df$group_sequence,"G")
dim(motif_df)
y=motif_df %>% 
  group_by(is_sd, quartile, dist) %>%
  summarise(mat=list(matrix(unlist(str_split(group_sequence,"")), ncol = 11, byrow = TRUE)))# %>%

z=y %>% 
 rowwise() %>% mutate(A = list(
   matrix(
     c(colSums(mat=="A"), 
     colSums(mat == "C"),
     colSums(mat == "G"),
     colSums(mat == "T")), 
     ncol=11, byrow=T, dimnames = list(c("A","C","G","T"),seq(1,11)) )
   )
   ) 

motif_df$m_len = sapply(motif_df$group_sequence,nchar)

library(ggseqlogo)
ggplot(data=z) + geom_logo(z$A)+theme_logo() + 
  facet_grid(is_sd*dist~quartile) 


#lm.df =  melt(data.table(gene.df), id.vars=c("gene","n_transcripts", "cut","is_sd", "n_cpg"), measure="methylated_frequency")
#lm.wide.df = data.table(lm.df %>% pivot_wider(names_from=cut, values_from=value, id_cols = c("gene","n_transcripts","is_sd")))
lm.wide.df = gene.df %>% 
  filter(min > -0.1 & max < 0.1) %>%
  pivot_wider(names_from=cut, values_from=methylated_frequency, id_cols = c("gene","n_transcripts","is_sd","genewidth"))
lm.wide.df
#bad = lm.wide.df[apply(lm.wide.df[,"(-1,-0.984]"],1,function(x)  length(unlist(x))) > 1]
#bad

# load package
library(sjPlot)
library(sjmisc)
library(sjlabelled)

m1 = lm(n_transcripts ~ . - gene -n_transcripts , data=lm.wide.df)
tab_model(m1)
}

