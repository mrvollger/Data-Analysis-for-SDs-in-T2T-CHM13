#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./plotutils.R")
# colors to use 
GRAY = "#2F4F4F"	
RED = "#af0404"
BLUE = "#3282b8"
NEWCOLOR = RED
OLDCOLOR = GRAY 
# output directories
FIGURE_DIR="~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Figures"
SUPP="~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Figures/Misc"
TABLES="~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Misc"
LOCAL_DATA="~/Desktop/Rdata"

if(F){
df = fread("../what_intersects_sv_flanks/annotated.sv.flanks.bed")
#df = fread("../what_intersects_sv_flanks/annotated.svs.bed")
freqs = fread("../what_intersects_sv_flanks/types.tbl") %>% filter(total_bp > 5e6) %>% 
  add_row(annotation="SD", total_bp=207556333) %>%
  mutate(annotation = factor(annotation, levels = unique(c("SD", sort(annotation)) )))
df = df %>% 
  filter(SVTYPE != ".") %>%
  mutate( annotation = SD ) %>%
  merge(freqs) %>%
  mutate(annotation = factor(annotation, levels = levels(freqs$annotation) ))  %>%
  filter(annotation != "SD") %>%
  data.table()

a = ggplot(data = df, aes(fill=SVTYPE)) +
  geom_bar(aes(x=annotation)) +
  facet_col(~SVTYPE, scales = "free_y") +
  theme_cowplot() +
  scale_fill_manual(values = c(NEWCOLOR,"black",BLUE)) +
  theme(axis.text.x = element_text(angle = -60, vjust=1, hjust=0),
        legend.position = "none") + 
  ggtitle("Annotations at SV flanks") + xlab("")

b = ggplot(data = df, aes(fill=SVTYPE)) +
  geom_bar(aes(x=annotation, weight = 1/total_bp)) +
  facet_col(~SVTYPE, scales = "free_y") +
  theme_cowplot() +
  scale_fill_manual(values = c(NEWCOLOR,"black",BLUE))+
  theme(axis.text.x = element_text(angle = -60, vjust=1, hjust=0),
        legend.position = "none") + 
  ggtitle("Relative enrichment of annotations at SV flanks")+ xlab("")

c=ggplot(data=freqs) + geom_bar(aes(annotation, weight=total_bp/1e6)) +   
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = -60, vjust=1, hjust=0),
        legend.position = "none") + 
  scale_y_continuous(label=comma) +
  ylab("Total Mbp in\n the genome")+ xlab("")

d=ggplot(data=df) + geom_bar(aes(annotation)) +   
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = -60, vjust=1, hjust=0),
        legend.position = "none") + 
  scale_y_continuous(label=comma) +
  ylab("Total counts in\n the flanks")+ xlab("")

#plot_grid(plot_grid(a,b, labels = "auto"),
#          c, ncol=1, labels = c(NA,"c"), rel_heights = c(2,1) )
fig = plot_grid(a,b,d,c, labels="auto", rel_heights = c(2,1))
ggsave(glue("{SUPP}/what_flanks_svs.pdf"), plot=fig, height=10, width=14)
}


make_sd_f <- function(pre, intra, flank = 100){
  #pre = "chm13_v1.0"
  total_sd_flank = fread(glue("../what_intersects_sv_flanks/{pre}.intra_{intra}.sd.{flank}.bed")) %>% summarise(sum(V3-V2))
  freqs = fread("../what_intersects_sv_flanks/types.sub.sub.tbl") %>% 
    #rbind(fread("../what_intersects_sv_flanks/types.sub.tbl"), fill=T) %>%
    #rbind(fread("../what_intersects_sv_flanks/types.tbl"), fill=T) %>%
    mutate(fraction=total_bp/3.11e9) %>%
    data.table()
  freqs
  
  sd_f = fread(glue("../what_intersects_sv_flanks/{pre}.intra_{intra}.sd.r.seq.{flank}.bed"))  %>% 
    merge(freqs) %>%
    data.table()
  
  sd_f$SD_track = pre

  enrichment_df = sd_f %>% group_by(SD_track, annotation, a2, a3) %>%
    summarise( observed = 100*sum(bp)/total_sd_flank[[1]], 
               expected = 100*unique(fraction),
               count = n(),
               seqs = list(seq),
               FlankGC = sum(str_count(seq, pattern = "G|g|C|c")) / sum(str_count(seq, pattern = "[^N]"))
               ) %>%
    ungroup() %>%
    #arrange(expected-observed) %>%
    mutate(CR=paste(annotation,a2,a3, sep=" / ")) %>%
    #arrange( - (observed - expected) ) %>% 
    #mutate(CR=factor(CR, levels = unique(CR)))
    mutate(enriched = log2(observed/expected) > 0 ) %>%
    mutate(status = case_when(intra=="1" ~ "intrachromosomal",
                              intra=="0" ~ "interchromosomal"
                              )
           ) %>%
    filter( (observed > 0.1) | (expected > 0.1) ) %>%
    #filter( (a2 %in% c("Alu")) |  (annotation %in% c("Satellite", "Low_complexity")) ) 
    data.table()
  enrichment_df
}

enrichment_df = rbind(make_sd_f("chm13_v1.0", "1"),
                      make_sd_f("chm13_v1.0", "0"),
                      make_sd_f("hg38_wgac", "1"),
                      make_sd_f("hg38_wgac", "0"),
                      #make_sd_f("hg38_wgac_merged","1"),
                      #make_sd_f("hg38_wgac_merged","0")
                      make_sd_f("hg38_sedef", "1"),
                      make_sd_f("hg38_sedef", "0")
                      )

sd_a = ggplot(data = enrichment_df %>%
                pivot_longer(cols = c("observed", "expected")),
                aes(fill=name)) +
  geom_col(aes(x=CR, y=value), position = "dodge") +
  theme_minimal_grid() +
  facet_grid(SD_track~status)+
  scale_fill_manual(values = c(OLDCOLOR, NEWCOLOR))+
  theme(axis.text.x = element_text(angle = -60, vjust=1, hjust=0),
        legend.position = "top", legend.title = element_blank()) + 
  ggtitle("Annotations at SD flanks")+ xlab("") +
  ylab("% of SD flank"); 

sd_b = ggplot(data = enrichment_df, 
              aes(y=log2(observed/expected),
                  fill=enriched)
              ) +
  geom_col(aes(x=CR), position = "dodge") +
  geom_text(aes(label = count, x=CR)) +
  theme_minimal_grid() +
  facet_grid(SD_track~status)+
  scale_fill_manual(values = c(OLDCOLOR, NEWCOLOR))+
  theme(axis.text.x = element_text(angle = -60, vjust=1, hjust=0),
        legend.position = "top") + 
  ggtitle("Log2 fold change at SD flanks")+ xlab("") +
  ylab("log2 fold change (observed/expected)"); 



sd_fig = plot_grid(sd_a, sd_b, ncol=2, align = "h")
#ggsave(glue("{SUPP}/what_flanks_sds_sat.pdf"), plot=sd_fig, height=12, width=24)
 



#
#
# find enirchment over hg38
#
#
en_13 = enrichment_df[enrichment_df$SD_track == "chm13_v1.0",]
en_38 = enrichment_df[enrichment_df$SD_track == "hg38_sedef",]
sd_c_df = merge(en_13, en_38, by=c("annotation", "a2", "a3", "CR", "status"), all.x = T)
toswap = is.na(sd_c_df$SD_track.y)
sd_c_df$SD_track.y[toswap] = "hg38_sedef" #sd_c_df$SD_track.x[toswap] 
sd_c_df[toswap, "count.y"] = 0
sd_c_df$logfold = log2(sd_c_df$count.x) - log2(sd_c_df$count.y)
sd_c_df$enriched = sd_c_df$logfold >= 0
sd_c_df$log_enriched = log2(sd_c_df$observed.x/sd_c_df$expected.x)  
sd_c_df$expected_count = sd_c_df$count.x * sd_c_df$expected.x/sd_c_df$observed.x

sd_c_df_filt = sd_c_df %>% filter( count.x >= 50, logfold >= 1, log_enriched > 1)

sd_c = ggplot(data = sd_c_df_filt, 
              aes(y=logfold)
              ) +
  geom_col(aes(x=CR, alpha = observed.x/expected.x),
           position = "dodge",
           fill=NEWCOLOR) +
  geom_text(aes(label = round(FlankGC.x*100), x=CR, y=-1.0, color="%GC")) +
  geom_text(aes(label = round(expected_count), x=CR, y=-0.75, color="Expected in T2T")) +
  geom_text(aes(label = count.y, x=CR, y=-0.5, color="GRCh38")) +
  geom_text(aes(label = count.x, x=CR, y=-0.25, color="T2T CHM13")) +
  facet_grid(status~.) +
  theme_minimal_grid() +
  scale_fill_continuous()+
  scale_alpha_continuous(trans="log2", range = c(0.1, 1))+
  scale_color_manual(values = c("darkgreen","black", OLDCOLOR, NEWCOLOR))+
  theme(axis.text.x = element_text(angle = -60, vjust=1, hjust=0),
        legend.position = "right") +
  guides(alpha=guide_legend("Fold enriched\nat SD flanks\ncompared to\nthe genome"),
         color=guide_legend("Raw SD\nflank counts")) +
  ylab("Log2 fold increase in annotations at\nSD flanks in T2T CHM13 vs GRCh38"); sd_c

ggsave(glue("{SUPP}/what_is_enriched_at_flanks.pdf"), 
       plot=sd_c, 
       height=8, width=length(unique(sd_c_df_filt$CR)), 
       limitsize = FALSE)


pdf(glue("{SUPP}/what_is_enriched_at_flanks_motifs.pdf"), height = 12, width = 12)
library(rGADEM); library("ggplotify"); library(ggseqlogo)
get_prime_motif<-function(row){
  seqs = row$seqs.x #sd_c_df_filt$seqs.x[2]
  dna = DNAStringSet(unlist(seqs))
  
  gadem<-GADEM(dna, nmotifs = 5)
  n = nOccurrences(gadem)
  pwm = rGADEM::getPWM(gadem)

  if(length(pwm) > 0){
    names(pwm) = paste(
      paste(n, length(dna), sep="/"),
      names(pwm),
      sep = "\n"
    )
    pwm = pwm[1:min(5, length(pwm))]
    ggplot() + 
      geom_logo(pwm) + 
      theme_logo() + 
      theme_cowplot() +
      ggtitle(row$CR) +
      facet_wrap(~seq_group, ncol=1, scales='free_x')
  } else {
    NULL
  }
}
apply(sd_c_df_filt, 1, get_prime_motif)
dev.off()
dev.off()


new_sd = SEDEF[overlaps(SEDEF, NEW)]# & overlaps(SEDEF[,10:12], NEW)]
bedlength(new_sd)/1e6; dim(new_sd)
fwrite(new_sd, file = "../what_intersects_sv_flanks/new_sds.bed", 
       row.names = F, col.names = F, sep="\t", scipen = 10000)











#
#
# read the meme resutls 
#
#
require(XML)
data <- xmlParse("~/Desktop/EichlerVolumes/mvollger/public_html/share/meme/simple_chm13_inter/meme.xml")
#data = xmlParse("https://eichlerlab.gs.washington.edu/help/mvollger/share/meme/simple_chm13_inter/meme.xml")
xml_data <- xmlToList(data)
  xml_data$scanned_sites_summary






