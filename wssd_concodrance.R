setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
theme_set(theme_cowplot())

if(F){
  non_rd = as.data.table(reduce(toGRanges(SEDEF)))
  kmer_cn = readbed("../wssd/v1.0/kmer/bed/chm13_wssd.sort.bed", tag="kmer CN")
  read_cn = readbed("../wssd/SD_Freeze1/bed/CHM13_wssd.sort.bed", tag="read CN")
  kmer_cn = kmer_cn[overlaps(kmer_cn, non_rd, mincov=0.99)]
  read_cn = read_cn[overlaps(read_cn, non_rd, mincov=0.99)]
  read_cn$CN_read = read_cn$V10
  kmer_cn$CN_kmer = kmer_cn$V10
  cn_set = merge(kmer_cn, read_cn, by=c("chr","start","end"))
  
  
  kmer_02 = readbed("../wssd/v1.0/kmer/bed/HG002_wssd.sort.bed", tag="kmer CN")
  read_02 = readbed("../wssd/SD_Freeze1/bed/HG002_wssd.sort.bed", tag="read CN")
  kmer_02 = kmer_02[overlaps(kmer_02, non_rd, mincov=0.99)]
  read_02 = read_02[overlaps(read_02, non_rd, mincov=0.99)]
  read_02$CN_read = read_02$V10
  kmer_02$CN_kmer = kmer_02$V10
  cn_02 = merge(read_02, kmer_02, by=c("chr","start","end"))
  
  
  kmer_733 = readbed("../wssd/v1.0/kmer/bed/HG00733_wssd.sort.bed", tag="kmer CN")
  read_733 = readbed("../wssd/SD_Freeze1/bed/HG00733_wssd.sort.bed", tag="read CN")
  kmer_733 = kmer_733[overlaps(kmer_733, non_rd, mincov=0.99)]
  read_733 = read_733[overlaps(read_733, non_rd, mincov=0.99)]
  read_733$CN_read = read_733$V10
  kmer_733$CN_kmer = kmer_733$V10
  cn_733 = merge(read_733, kmer_733, by=c("chr","start","end"))
}
dim(read_cn)
dim(kmer_cn)




ggplot(data = cn_set, aes(x=(CN_kmer+1), y = (CN_read + 1)) 
       ) +
  geom_hex(bins = 300)+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")


cor

cor.test(cn_set$CN_read, cn_set$CN_kmer, method=c("pearson"))
cor.test(cn_02$CN_read, cn_02$CN_kmer, method=c("pearson"))
