#!/usr/bin/env Rscript
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")
library(ggridges)

peri = readbed("../random_stat_scripts/split_sds/sds.peri.dedup.bed", "pericentromeric")
intra = readbed("../random_stat_scripts/split_sds/sds.intra.dedup.bed", "intrachromosomal")
inter = readbed("../random_stat_scripts/split_sds/sds.inter.dedup.bed", "interchromosomal")
acro = readbed("../random_stat_scripts/split_sds/sds.acro.dedup.bed", "acrocentric")
telo = readbed("../random_stat_scripts/split_sds/sds.telo.dedup.bed", "subtelomeric")

df = rbind(intra, inter, acro, peri, telo)
x = "#chr1   start1  end1    name    fakeScore       strand1 start1  end1    color   chr2    start2  end2    score   strand2 max_len aln_len indel_a indel_b alnB    matchB  mismatchB transitionsB    transversions   fracMatch       fracMatchIndel  jck     k2K     aln_gaps        uppercaseA      uppercaseB      uppercaseMatches        aln_matches       aln_mismatches  aln_gap_bases   filter_score    sat_bases"
header = strsplit(x, "\\s+")[[1]]
colnames(df)[4:length(header)] = header[4:length(header)]

names = rev(sort(unique(df$Assembly)))
df$Assembly = factor(df$Assembly, levels = names )
colors = c("green", OLDCOLOR, BLUE, NEWCOLOR, "#FF0000")
names(colors) = names

ggplot(data=df) +
  geom_point(aes(x=max_len, y = fracMatch), alpha=0.8, size=0.05) +
  geom_density2d(aes(x=max_len, y = fracMatch), size=1, alpha=0.8) +
  #scale_x_continuous(trans = "log10", labels = comma) +
  #scale_fill_brewer(palette ="Spectral", direction = -1)+
  facet_wrap(vars(Assembly)) +
  theme_cowplot()



p1 = ggplot(data=df)+
  geom_density_ridges(aes(max_len, y=Assembly, fill=Assembly))+
  xlab("Segmental duplication length (bp)")+
  scale_fill_manual( values = colors) +
  theme_cowplot()+
  theme(legend.position = "none") + 
  scale_x_continuous(trans = "log10", labels = comma) #facet_zoom( x = max_len >= 5e4 & max_len <= 3e5, zoom.size=2, horizontal = F) ;p1

p2 = ggplot(data=df)+
  geom_density_ridges(aes(fracMatch*100, y=Assembly, fill=Assembly))+
  xlab("Percent Identity")+
  scale_x_continuous(breaks = seq(90,100))+
  scale_fill_manual( values = colors) +
  theme_cowplot()+
  theme(legend.position = "none")

p = plot_grid(p1,p2)

ggsave("supp/sd_length_identity.pdf", plot = p, width = 16, height = 8)
p

fileConn<-file("supp/sd_length_identity.txt")
out = c()
for(n in names){
  for(j in names){
    if(n!=j){
      result =  wilcox.test(df[Assembly==n]$max_len, df[Assembly==j]$max_len, alternative = "greater")
      if(result$p.value<0.5){
        #print(c(n,">",j))
        #print(result$p.value)
        out = c(out, c(paste(n,"<",j), result$p.value))
        }  
      }
  }
}
writeLines(out, fileConn)
close(fileConn)
