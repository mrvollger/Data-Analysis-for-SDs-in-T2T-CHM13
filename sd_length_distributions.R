#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")
library(ggridges)

tmp = copy(SEDEF[SEDEF$original])
peri = tmp[peri>0 | peri2 > 0]; peri$Assembly = "pericentromeric"
acro = tmp[acro>0 | acro2 > 0]; acro$Assembly = "acrocentric"
telo = tmp[telo>0 | telo2 > 0]; telo$Assembly = "telomeric"
intra = tmp[chr == chr2]; intra$Assembly = "intrachromosomal"
inter = tmp[chr != chr2]; inter$Assembly = "interchromosomal"
notacro = tmp[acro == 0 & acro2 == 0]; notacro$Assembly = "not acrocentric"
all = tmp; all$Assembly = "All"
df = rbind(intra, inter, acro, peri, telo, all, notacro)
df$length = df$matchB + df$mismatchB #df$end - df$start
df$Assembly = factor(df$Assembly, levels =  unique(df$Assembly))
df$Assembly

x = "#chr1   start1  end1    name    fakeScore       strand1 start1  end1    color   chr2    start2  end2    score   strand2 max_len aln_len indel_a indel_b alnB    matchB  mismatchB transitionsB    transversions   fracMatch       fracMatchIndel  jck     k2K     aln_gaps        uppercaseA      uppercaseB      uppercaseMatches        aln_matches       aln_mismatches  aln_gap_bases   filter_score    sat_bases"
header = strsplit(x, "\\s+")[[1]]
colnames(df)[4:length(header)] = header[4:length(header)]

names = levels(df$Assembly)
colors = c(OLDCOLOR, BLUE, NEWCOLOR, "#FF0000","green", "black", 'blue')
names(colors) = names

ggplot(data=df) +
  geom_point(aes(x=max_len, y = fracMatch), alpha=0.8, size=0.05) +
  geom_density2d(aes(x=max_len, y = fracMatch), size=1, alpha=0.8) +
  #scale_x_continuous(trans = "log10", labels = comma) +
  #scale_fill_brewer(palette ="Spectral", direction = -1)+
  facet_wrap(vars(Assembly)) +
  theme_cowplot()


info.df = df %>% group_by(Assembly) %>% 
  group_modify(~length_stats(.x)) %>% 
  mutate(label=as.factor(label))
info.df$y = rep(seq(.1,.8,.7/(length(unique(info.df$label))-1)), length(unique(info.df$Assembly)))
info.df

p1 = ggplot(data=df)+
  geom_density_ridges(aes(max_len, y=Assembly, fill=Assembly))+
  xlab("Segmental duplication length (bp)")+
  geom_text(data=info.df, 
            aes(label=paste(label, round(value)), x= 1e6, y=as.numeric(Assembly)+y, group = Assembly, color=Assembly),
            #direction = "y",
            #ylim = c(0,NA),
            hjust = 0)+
  scale_fill_manual( values = colors) +
  scale_color_manual( values = colors) +
  theme_cowplot()+
  theme(legend.position = "none") + 
  scale_x_continuous(trans = "log10", labels = comma) #facet_zoom( x = max_len >= 5e4 & max_len <= 3e5, zoom.size=2, horizontal = F) ;p1
p1
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

tmp

length_stats(df[Assembly=="acrocentric"])
length_stats()
