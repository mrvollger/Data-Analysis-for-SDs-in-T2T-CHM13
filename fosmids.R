#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#source("plotutils.R")


tested_sds = readbed("../FES_stuff/acro_cn/probe.mappings.bed", "sds")
probe = c("174222_ABC10_2_1_000044559800_C15","174222_ABC10_2_1_000044587300_G6","171417_ABC10_2_1_000045520200_K20")
#tested_sds=tested_sds[probe_id %in% probe]
n =length(unique(tested_sds$probe_id))



#
# read in fes results
#
results = fread("../FES_stuff/acro_cn/results.csv")
results = melt(results, id.vars =c("sample", "CLONES"))
results$loc = gsub("\\(het\\)", "", results$value)
results$chr = results$variable
results= results[!is.na(results$loc)]
results=results[loc !="" ]
results

cens = CYTO[CYTO$gieStain == "gpos100"]
pcen = resize(cens, 5000000, fix="end")
qcen = resize(cens, 5000000, fix="start")

a = DataFrame(chr = as.character(seqnames(pcen)), start = start(ranges(pcen)), end = end(ranges(pcen)) ); a$loc = "pcen"
b = DataFrame(chr = as.character(seqnames(qcen)), start = start(ranges(qcen)), end = end(ranges(qcen)) ); b$loc = "qcen"
c = DataFrame(chr = FAI$chr, start = FAI$chrlen-5000000, end = FAI$chrlen); c$loc = "qter"
d = DataFrame(chr = as.character(seqnames(qcen)), start = 0, end = 5000000); d$loc = "pter"
e = DataFrame(chr=("chr4"), start = c(75300001), end = c(87100000), loc = "q21")
f = DataFrame(chr = as.character(seqnames(qcen)), start = start(ranges(qcen)), end = end(ranges(qcen)) ); f$loc = "qpcen"
convert_results = as.data.table(rbind(a,b,c,d,e,f))


fish_results = merge(results, convert_results, by = c("chr", "loc"))
fish_results= fish_results[chr %in% NOM]
fish_results <- as.data.table(merge(as.data.frame(fish_results), 
                                   data.table(color=brewer.pal(length(unique(fish_results$sample)), "Set3"), 
                                              sample = unique(fish_results$sample))))
gr_fish = toGRanges(data.table(fish_results$chr, fish_results$start, fish_results$end))

#
# show only the tested probes
#
all = tested_sds[tested_sds$probe_id %in% fish_results$CLONES]
probe = unique(all$probe_id)
n=length(probe)
#
# plot results
#

myplots <- vector('list', n)

for(i in seq(n)){
  df <- all[probe_id==probe[i]]
  df2 <- fish_results[CLONES==probe[i]]
  chromosomes <- CHRS[ CHRS %in% df$chr | CHRS %in% df2$chr ]
  #df <- df[df$end-df$start > 9999]
  
  myplots[[i]] = as.ggplot(expression(
    kp <- plotKaryotype(genome = GENOME, cytobands = CYTO, chromosomes = chromosomes, plot.type = 2),
    kpAddMainTitle(kp, main=probe[i]),
    kpPlotRegions(kp, data = reduce(toGRanges(data.table(chr = as.character(df$chr), start = df$start, end=df$end))+1000000), 
                  col = NEWCOLOR),
    kpPlotRegions(kp, data = toGRanges(data.table(chr = as.character(df2$chr), start = df2$start, end=df2$end))+5000000, 
                  col = df2$color, border = df2$color, data.panel = 2)
  ))
}
myplots[[1]]


pdf("~/Desktop/fosmid.pdf", height = 4, width = 16)
for(i in seq(n)){
  print(myplots[[i]])
}
dev.off()




