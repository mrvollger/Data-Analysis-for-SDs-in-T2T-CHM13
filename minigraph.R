#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if(! require("tidyverse")) install.packages("tidyverse")
if(! require("ggnewscale")) install.packages("ggnewscale")
if(! require("ggrepel")) install.packages("ggrepel")
if(! require("data.table")) install.packages("data.table")
if(! require("glue")) install.packages("glue")
if(! require("RColorBrewer")) install.packages("RColorBrewer")
if(! require("scales")) install.packages("scales")
if(! require("cowplot")) install.packages("cowplot")
source("minigraph_utils.R")


make_minigraph_figure = function(gene, search_gene, simple_gene){
  #
  # load data, and make names human readable
  # 
  DF = fread(glue("{DATA_DIR}/Minigraph/{gene}.tbl"))
  DUPLICONS = tri_bed(glue("{DATA_DIR}/Masked/{gene}_dupmasker_colors.bed"))
  All_genes = fread(glue("{DATA_DIR}/Liftoff/{gene}.all.bed"))
  names(All_genes)[1:3] = c("chr","start","end")
  GFA = glue("{DATA_DIR}/Minigraph/{gene}.gfa")
  CSV = fread(glue("{DATA_DIR}/Minigraph/{gene}.csv"))
  
  #
  # Filter data to only show top results
  # 
  R_NAMES = unique(DF$r)
  Q_NAMES = unique(DF$q)
  Q_ONLY = Q_NAMES[!Q_NAMES %in% R_NAMES]
  NUM_R = length(R_NAMES)
  MAX_SM = max(NUM_R, MAX_ROWS)
  KEEP = as.factor(c(as.character(R_NAMES), as.character(Q_ONLY[1:(MAX_SM-NUM_R)])))
  KEEP = KEEP[!is.na(KEEP)]
  
  CHM13_name = grep("CHM13.pri", c(R_NAMES, Q_NAMES), value = T)[1]
  GRCh38_name = grep("GRCh38chrOnly", c(R_NAMES, Q_NAMES), value = T)[1]
  
  # orgganize by length
  KEEP = unique(c(CHM13_name, GRCh38_name, merge(data.table(q=KEEP), DF[, c("q","ql")])[order(ql)]$q))
  
  # cleanup the DF
  df = clean_df(DF, KEEP)
  duplicons = clean_df(DUPLICONS, KEEP)
  all_genes = clean_df(All_genes, KEEP)
  csv = clean_df(CSV, KEEP)
  
  #
  # code
  #
  
  df$y = as.numeric(df$q)
  df$y = 1
  df$id = 1:nrow(df)
  tmp = data.table(df[,q,ql] %>% group_by(q, ql) %>% unique())
  names(tmp)  = c("rl","r")
  df = merge(df, tmp, by="r")
  
  
  
  #
  # make a duplicons legend
  #
  all_genes$q = all_genes$chr
  all_genes$target  = grepl(gsub("_.*","", search_gene), all_genes$V4) 
  all_genes$arrow = "last"
  all_genes[V6=="-"]$arrow = "first"
  all_genes = all_genes[chr %in% unique(df$q)]
  all_genes$length = all_genes$end - all_genes$start
  # filter to only largest per, and then filter by size
  all_genes = all_genes[all_genes[, .I[which.max(length)], by=c("chr","V4","V6")]$V1]
  #all_genes = all_genes[length>1000]
  
  genes = all_genes[all_genes$target]
  #gene_rgn = genes[which.max(genes$end-genes$start)][1]
  if(F){
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
  }else{
    dup_legend = ggplot()
  }
  
  
  #
  #  set up the colors
  #
  otherr =  unique(df$r)[ !grepl("CHM13|GRCh38", unique(df$r), ignore.case = TRUE, perl = TRUE)]
  if( length(unique(df$r))==1 ){
    rcolors =  NEWCOLOR
  }else if( length(unique(df$r))==2 ){
    rcolors=c(NEWCOLOR, OLDCOLOR)
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
  
  #
  # set up the dataframe for ploting
  #
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
  
  # add some labels for the figure
  label.df = data.frame(r = plot.df$r[1], q = plot.df$r[1],
                        y = c(1+dup_offset/2, 1-(dup_offset-s)/2, 1-1*dup_offset),
                        x = 1.02 * max(plot.df$ql[1]),
                        label = c("Squashed dot plot", "Genes", "Duplicons")
                        )
  
  #
  # make the figure
  #
  
  # establish the x and y axis, add diagnonal dot plot
  p1 <- ggplot() + 
    geom_segment(data=lines, aes(y=y, yend=y, x=x, xend=xend, color=q), size=1, alpha=0.9)+
    geom_segment(data = lines, aes(x=-max(xend)/50, xend = -max(xend)/50, y=0-s/2, yend=1+s, color=q), size=4)+
    geom_polygon(data=plot.df, aes(x=xs, y=ys, group=id, fill=r), color="black", alpha=0.9)+
    geom_text(data = label.df, aes(x=x, y=y, label = label, group=r), hjust = 0, vjust = 0.5, fontface =2)+
    coord_cartesian(clip="off") +
    theme_cowplot()
    
  #+ # turn off clipping for text labels
    #ggtitle(glue("Variation in {simple_gene}"))
  
  bases_of_haplotypes = sum(unique(df[!q %in% c("GRCh38"), q, ql])$ql)
  # add the text for SV events
  tmpdf = df %>% filter(!grepl("_0kbp|Syntenic", description)) %>% mutate(description = tolower(str_replace(description, "_", " ")))
  labdf=tmpdf %>% 
    separate(description, c("type", "size"), sep=" ") %>%
    mutate(size = as.numeric(str_replace(size,"kbp", ""))) %>% 
    group_by(q, type) %>%
    summarise(count = n(), kbp=sum(size), ql=max(ql)) %>%
    mutate(tmp = paste(paste(type, kbp, sep=": "), glue("({count})"), sep=" kbp ")   ) %>%
    ungroup() %>% group_by(q) %>%
    summarise(label = str_c(tmp, collapse = '\n'), ql=max(ql))
      
  tmpdf %>% mutate( SV_kbp = as.numeric(gsub("ins|inv|del|kbp| ", "", description)),
                    SV_type = gsub(" .+", "", description)) %>% 
    #filter(q != "GRCh38") %>%
    group_by(SV_type) %>%
    summarise(kbp = sum(SV_kbp), frac = sum(SV_kbp)*1000/bases_of_haplotypes*100 )
  
  #
  # fancy labels or summary table
  #
  if(T){
    p1 = p1 + geom_label(data=labdf,
                               aes(x=max(ql), y=1, label=label), 
                         hjust=0, vjust=0.5)
  } else if(nrow(tmpdf)>0){
    p1 = p1 + geom_label_repel(data=tmpdf, 
                    aes(x=(qe+qs)/2, y=1, label = description),#, fill=r),
                    #fontface = "bold",
                    ylim         = c(NA, .9),
                    direction    = "x",
                    angle        = 0,
                    segment.size = 0.5 ) 
  } 
  
  # continue ploting
  p = p1 +
    # clear previous scales
    scale_y_continuous(limits=c(NA,NA))+
    scale_fill_manual(values=colors)+ 
    scale_color_manual(values=colors) +
    new_scale_fill() + new_scale_color() +
    
    # add in the gene data
    geom_segment(data=all_genes %>% filter(arrow == "first"), 
             aes(x=start, xend=end,
                 y=1 - (dup_offset-s*target)/2, yend=1 - (dup_offset-s*target)/2,
                 color=target, group=chr),
             arrow = arrow(length = unit(0.08, "npc"), ends = "first"),
             size=.5)+
    geom_segment(data=all_genes %>% filter(arrow == "last"), 
                 aes(x=start, xend=end,
                     y=1 - (dup_offset-s*target)/2, yend=1 - (dup_offset-s*target)/2,
                     color=target, group=chr),
                 arrow = arrow(length = unit(0.08, "npc"), ends = "last"),
                 size=.5)+
    scale_color_manual(values = c(`FALSE`="black", `TRUE`="red")) +
    new_scale_color() + 
    
    # add in the duplicons
    geom_polygon(data=duplicons, aes(
      x=xs, y=ys-dup_offset, group=tri_id, fill=color
    )) + 
    scale_fill_manual(values=dcolor )+
    
    # add figure astetics
    ylab("")+#ylab("Haplotype resolved assemblies") + 
    xlab("Genomic position (bp)")+
    #scale_y_continuous( labels = levels(df$q), breaks=1:length(levels(df$q)) ) + 
    scale_x_continuous(labels = comma, expand = c(0,NA))+
    facet_wrap(vars(q), strip.position = "left", ncol=1) +
    theme_cowplot() +
    theme(legend.position = "none",
          plot.title = element_text(hjust=0.5),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
          strip.background = element_rect(colour=NA, fill=NA),
          panel.border = element_rect(color = rgb(col2rgb("gray")[,1]["red"],
                                                  col2rgb("gray")[,1]["green"],
                                                  col2rgb("gray")[,1]["blue"],
                                                  0.5, maxColorValue = 255), fill = NA, size = .1),
          plot.margin = unit(c(1,3,1,1), "cm"),
          strip.text.y.left   = element_text(angle=0, hjust=0)
          ); p
  
  
    
  #
  #
  # Make the bandage figure
  #
  #
  if( ! dir.exists("temp")){
    dir.create("temp")
  }
  colorcsv = merge(csv, data.table(Label = names(colors), Color=colors), by="Label")[,c("Name","Color")]
  write.table(colorcsv, file="temp/tmp.gfa.csv", sep=",", quote = FALSE, row.names = FALSE)
  
  # ./Bandage.app/Contents/MacOS/Bandage image ~/Desktop/EichlerVolumes/chm13_t2t/nobackups/minigraph/loci_of_interest/CHR1_QCEN.gfa ~/Desktop/tmp.png --colors ~/Desktop/EichlerVolumes/chm13_t2t/nobackups/plots/temp/tmp.gfa.csv
  gfa_h = 400*length(unique(df$q))
  
  system(glue("{BANDAGE_PATH}/Bandage image {GFA} temp/{gene}.png --linear --height {gfa_h}  --colors temp/tmp.gfa.csv"))
  system(glue("convert temp/{gene}.png -rotate 90 temp/{gene}.r.png"))
  #system(glue("rm temp/{gene}.png temp/tmp.gfa.csv"))
  graph_figure = ggdraw() + draw_image(glue("temp/{gene}.r.png"))
  
  
  #
  #
  # make a legned for the figure
  #
  #
  num_q_seqs = length(unique(df$q))
  
  for_l = ggplot(data=df) + 
    geom_bar(aes(x=q, fill=q)) + scale_fill_manual(values = colors)+ theme_cowplot()+
    theme(legend.position = "bottom", 
          legend.justification = "center",
          plot.title = element_text(hjust=0.5))+
    guides(fill=guide_legend(ncol=round(num_q_seqs),
                             title = "Assembly")
           ) +
    ggtitle(glue("Variation in {simple_gene}"))
  legend = cowplot::get_legend(for_l)
  title = cowplot::get_title(for_l)
  #plot_grid(title,legend)
  
  #
  #
  # make the final figure 
  #
  #
  if( ! dir.exists("minigraph_figures")){
    dir.create("minigraph_figures")
  }
  
  #figure=plot_grid(graph_figure, p, nrow=2, rel_heights = c(1, 1.5*length(unique(df$q)))); figure
  figure1=cowplot::plot_grid(graph_figure, p, ncol=2, rel_widths = c(2,4), labels = "auto")
  figure = cowplot::plot_grid(title, figure1, nrow=2, rel_heights = c(.5,max(num_q_seqs,8)))
  adjust = .8
  ggsave(glue("{FIGURE_DIR}/Minigraph Figures/{gene}.pdf"), height = adjust*1.25*max(num_q_seqs,8), width = adjust*20, plot = figure)

}




#
#  START READING
#

# change this if you want
MAX_ROWS = 150 # the maximum number of rows to plot, though it always plots all the seqs the contribute to building the graph
# CHANGE THIS !!!!!!!!!!!
DATA_DIR = "../orthology_analysis/TBC1D3/minigraph/pull_sd_regions"
DATA_DIR = "../sd_regions_in_hifi_wga/pull_by_regions_snake_results/pull_sd_regions"
BANDAGE_PATH = "~/software/Bandage//Bandage.app/Contents/MacOS/"
FIGURE_DIR="~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Figures/"

# READ HERE!!!!!!!!!!!!
#
# set up a gene name 
#
# gene = "region you name you used in pull_region"
# search_gene = "gene\pattern to search for in the gene file"
# simple_gene = "title for plot: variation in {simple_gene}"
make_minigraph_figure("TBC1D3_1", "TBC1D3", "TBC1D3 (1)")
make_minigraph_figure("TBC1D3_2", "TBC1D3", "TBC1D3 (2)")
#make_minigraph_figure("middle_TBC1D3_unique_sequence", "TBC1D3","between TBC1D3 site 1 and 2")
make_minigraph_figure("SMN", "SMN","SMN")
make_minigraph_figure("CYP2D6", "CYP2D","CYP2D6")
make_minigraph_figure("ARHGAP11", "ARHGAP11","ARHGAP11")
make_minigraph_figure("BOLA2_NPIPB", "BOLA2","BOLA2")
make_minigraph_figure("NOTCH2", "NOTCH2","NOTCH2")
make_minigraph_figure("NOTCH2NL", "NOTCH2NL","NOTCH2NL")
make_minigraph_figure("TCAF2", "TCAF","TCAF")
make_minigraph_figure("LPA", "LPA","LPA")
make_minigraph_figure("SRGAP2C","SRGAP2","SRGAP2")

# SRGAP2A has zero heterozygosity, ignore. 
#gene="SRGAP2A"; search_gene = "SRGAP2"; simple_gene="SRGAP2"
#make_minigraph_figure("SRGAP2A","SRGAP2","SRGAP2")
