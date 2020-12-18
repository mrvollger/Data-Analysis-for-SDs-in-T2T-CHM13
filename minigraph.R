#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#source("plotutils.R")
if(! require("tidyverse")) install.packages("tidyverse")
#library(tidyverse)
if(! require("ggnewscale")) install.packages("ggnewscale")
if(! require("ggrepel")) install.packages("ggrepel")
if(! require("data.table")) install.packages("data.table")
if(! require("glue")) install.packages("glue")
if(! require("RColorBrewer")) install.packages("RColorBrewer")
if(! require("scales")) install.packages("scales")
if(! require("cowplot")) install.packages("cowplot")

GRAY = "#2F4F4F"	
RED = "#af0404"
BLUE = "#3282b8"
NEWCOLOR = RED
OLDCOLOR = GRAY 


gene="SMN_LL"
gene="LPA"
gene="TCAF2"
gene="NOTCH2NLA"
gene="EIF3C"
gene="TBC1D3_1"
gene="CYP2D6"
gene="ARHGAP11_L"
gene="SRGAP2B_D"
gene="CHR1_QCEN_L"
gene="DEF_LLL"
gene="TBC1D3_2"

tri_bed <- function(f, s=.2, allowed_names=c()){
  #df = fread(glue("../sd_regions_in_hifi_wga/lpa/minimiro/temp_minimiro/{gene}_query.fasta.duplicons.extra")); df
  df=fread(f)
  names(df)[1:3]=c("chr","start","end")
  if(length(allowed_names) > 0){
    df = df[chr %in% allowed_names]
  }
  df = df[order(chr,start)]
  df$chr = factor(df$chr)
  df$y = as.numeric(df$chr)
  df$y = 1
  df$tri_id = 1:nrow(df)
  if("orient" %in% names(df)){
    df$strand = "+"
    df$strand[df$orient == "R"] = "-"
  }
  # allow to be faceted together
  df$q = df$chr
  zs = s/3
  data.table(df %>% 
    rowwise() %>% 
    mutate(xs=list(c(start, start, end, end)), 
           ys = case_when(
      strand == "+" ~ list(c(y+s,y-s,y-zs,y+zs)),
      strand == "-" ~ list(c(y-zs,y+zs,y+s,y-s))
    )) %>% unnest(cols=c("xs","ys")))
}

df = fread(glue("../minigraph/cyp2d6/{gene}.tbl"))
df = df[order(q,qs)]
df$q = factor(df$q)
df$r = factor(df$r)
df$y = as.numeric(df$q)
df$y = 1
df$id = 1:nrow(df)
tmp = data.table(df[,q,ql] %>% group_by(q, ql) %>% unique())
names(tmp)  = c("rl","r")
df = merge(df, tmp, by="r")

duplicons = tri_bed(glue("../sd_regions_in_hifi_wga/lpa/minimiro/temp_minimiro/{gene}_query.fasta.duplicons.extra"), allowed_names = unique(df$q, df$r))


#
# make a duplicons legend
#
genes = readbed(glue("../sd_regions_in_hifi_wga/lpa/minimiro/Liftoff/{gene}.all.bed"),"genes"); genes = genes[grep(gsub("_.*","",gene),genes$V4)] 
#gene_rgn = genes[which.max(genes$end-genes$start)][1]
gene_rgn=genes[genes$chr=="CHM13.pri__1"][1]
has_both = findOverlapPairs(toGRanges(duplicons) ,toGRanges(gene_rgn)+1000)
onedup = as.data.table(first(has_both))
mlen = max(df$qe)
dup_legend = ggplot()+
  geom_polygon(data=onedup, aes(x=xs-min(xs)+mlen/2, y=ys-.2, group=tri_id, fill=color))+
  geom_line( aes(x=c(0, max(df$qe)), y=c(1,1) ) , color="white")+
  scale_fill_identity()+
  theme_map()+
  labs(title=gsub("_.*","",gene))+theme(plot.title = element_text(hjust=0.5)) ;dup_legend
#
# 
#



#duplicons = duplicons[chr %in% unique(df$q)]

otherr =  unique(df$r)[ !grepl("CHM13|GRCh38", unique(df$r), ignore.case = TRUE, perl = TRUE)]
if( length(unique(df$r))==1 ){
  rcolors =  NEWCOLOR
}else if(length(unique(df$r))<13) {
  rcolors = brewer.pal(length(unique(df$r)), "Set3") 
} else{
  rcolors =  colorRampPalette(brewer.pal(12, "Set3"))(length(unique(df$r)))
}


names(rcolors) = unique(df$r)#[ grepl("CHM13|GRCh38", unique(df$r), ignore.case = TRUE, perl = TRUE)]
qnames = unique(df$q[! df$q %in% df$r ]) 
qcolors = rep("darkgray",length(qnames))
names(qcolors) = qnames
colors = c(rcolors,  qcolors)

colors[grepl("CHM13" , names(colors))] = NEWCOLOR
colors[grepl("GRCh38" , names(colors))] = OLDCOLOR

colors 

s=.6
dup_offset = 1

plot.df = df %>% 
  rowwise() %>% 
  mutate(xs=list(c(qs, qs, qe, qe)), ys = case_when(
      (strand == "+" & direction == ">") | (strand == "-" & direction == "<") ~ list(c(y , y + s*rs/rl, y + s*re/rl, y )),
      TRUE ~ list(c(y , y + s*re/rl, y + s*rs/rl, y ))
    )) %>% unnest(cols=c("xs","ys"))

lines = df %>% group_by(q) %>% summarise(y = (unique(y)-dup_offset), x=min(qs), xend=max(qe))

dcolor = unique(duplicons$color)
names(dcolor) = unique(duplicons$color)

p <- ggplot() + 

  geom_segment(data=lines, aes(y=y, yend=y, x=x, xend=xend, color=q), size=1, alpha=0.9)+
  geom_segment(data = lines, aes(x=-max(xend)/50, xend = -max(xend)/50, y=0-s/2, yend=1+s, color=q), size=4)

tmpdf = df %>% filter(!grepl("_0kbp|Syntenic", description)) %>% mutate(description = tolower(str_replace(description, "_", " ")))
if(nrow(tmpdf)>0){
  p = p + geom_label_repel(data=tmpdf, 
                  aes(x=(qe+qs)/2, y=1, label = description),#, fill=r),
                  #fontface = "bold",
                  ylim         = c(NA, .9),
                  direction    = "x",
                  angle        = 0,
                  segment.size = 0.5 ) +

  geom_polygon(data=plot.df, aes(x=xs, y=ys, group=id, fill=r), color="black", alpha=0.9)
}

p = p +
  scale_y_continuous(limits=c(NA,NA))+
  scale_fill_manual(values=colors)+ 
  scale_color_manual(values=colors) +
  new_scale_fill() + new_scale_color() + 
  
  geom_polygon(data=duplicons, aes(
    x=xs, y=ys-dup_offset, group=tri_id, fill=color
  )) + 
  scale_fill_manual(values=dcolor )+
    
  ylab("Haplotype resolved assemblies") + 
  xlab("genomic position (bp)")+
  #scale_y_continuous( labels = levels(df$q), breaks=1:length(levels(df$q)) ) + 
  scale_x_continuous(labels = comma)+
  facet_wrap(vars(q), nrow=length(unique(df$q))) +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        strip.background = element_rect(colour=NA, fill=NA)); p

if( ! dir.exists("minigraph_figures")){
  dir.create("minigraph_figures")
}

figure=plot_grid(dup_legend, p, nrow=2, rel_heights = c(1, 1.5*length(unique(df$q))))
ggsave(glue("minigraph_figures/{gene}.pdf"), height = 1.5*length(unique(df$q)), width = 16, plot = figure)
  

#cols <- function(a) image(1:length(a), 1, as.matrix(1:length(a)), col=a, axes=FALSE , xlab="", ylab="")
#cols(unique(dcolor))


