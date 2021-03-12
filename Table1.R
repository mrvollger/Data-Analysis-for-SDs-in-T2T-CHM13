#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(webshot)
library(openxlsx)
source("plotutils.R")
if( ! dir.exists("Tables")){
  dir.create("Tables")
}

hs <- createStyle(
  textDecoration = "BOLD", fontColour = "#000000", fontSize = 12,
  border="bottom", borderStyle="medium",halign="center"
  #fontName = "Arial Narrow", fgFill = "#4F80BD"
)

nonr <- function(df){
  grange = toGRanges(data.table(df))
  gr = GenomicRanges::reduce(grange)
  sum(width(gr))/1e6
}

my_sum <- function(df){
  nbp = nonr(df)
  #
  inter = nonr( df[df$chr!=df$chr2,] )
  ninter = length(unique( df[df$chr!=df$chr2,]$unique_id ))
  #
  intra = nonr( df[df$chr==df$chr2,] )
  nintra = length(unique( df[df$chr==df$chr2,]$unique_id ))
  #
  acro = nonr(df[df$acro>0,])
  nacro = length(unique( df[df$acro>0,]$unique_id ))
  #
  peri = nonr(df[df$peri>0,])
  nperi = length(unique( df[df$peri>0,]$unique_id ))
  #
  telo = nonr(df[df$telo>0,])
  ntelo = length(unique( df[df$telo>0,]$unique_id ))
  
  
  df %>% summarise( `Gbp` = sum(unique(chrlen))/1e9,
                    `% SD` = nbp/(sum(unique(chrlen))/1e6)*100,
                    `SD (Mbp)` = nbp,
                    `# SDs` = length(unique(unique_id))) %>%
    mutate(`inter (Mbp)` = inter, `# inter` = ninter,
           `intra (Mbp)` = intra, `# intra` = nintra,
           `acro (Mbp)` = acro, `# acro` = nacro,
           `peri (Mbp)` = peri , `# peri` = nperi,
           `telo (Mbp)` = telo, `# telo` = ntelo)
}

new_or_sv = copy(SEDEF[overlaps(SEDEF, NEW)])
new_or_sv$Assembly = "New or structurally variable"



out = rbind(SEDEF, SEDEF_38, new_or_sv) %>% 
  group_by(Assembly) %>%
  group_modify( ~ my_sum(.x)) %>%
  ungroup() %>% 
  arrange(-`SD (Mbp)`)

out[3,c("Gbp","% SD")] = list(nonr(NEW)/1e3, nonr(new_or_sv)/(nonr(NEW))*100)
diff = tibble(Assembly="Difference", out[1,2:length(out)] - out[2,2:length(out)] )
out=rbind(out,diff)


z="
HiFi coverage CN estimate
chr	2N	1N
13	150	75
14	39	19.5
15	52	44 (8)
21	124	62
22	44	22
ILMN k-mer CN estimate 
chr	2N	1N
13	168	84
14	38	19
15	50	42 (8)
21	108	54
22	36	18
"
unit=45e3
n = ceiling(75+19.5+44+62+22)
intra_n = ceiling(75*74/2 + 19.5*18.5/2 + 44*43/2 + 62*61/2 + 22*21/2)
rDNA = copy(out[1,])

# lead
rDNA$Assembly = "T2T CHM13 + rDNA (estimate)"
rDNA$Gbp =   0
rDNA$`SD (Mbp)` = unit * n / 1e6
rDNA$`# SDs` = n * ( n - 1 ) / 2
rDNA$`% SD` = rDNA$`SD (Mbp)`/( out$Gbp[1] * 1e3) *100


# bp
rDNA$`intra (Mbp)` = rDNA$`SD (Mbp)`
rDNA$`inter (Mbp)` = rDNA$`SD (Mbp)`
rDNA$`acro (Mbp)` = rDNA$`SD (Mbp)`
rDNA$`peri (Mbp)` = rDNA$`SD (Mbp)`
rDNA$`telo (Mbp)` = 0

# counts 
rDNA$`# intra` = intra_n
rDNA$`# inter` = rDNA$`# SDs` - rDNA$`# intra`
rDNA$`# acro` = rDNA$`# SDs`
rDNA$`# peri` = rDNA$`# SDs`
rDNA$`# telo` = 0

out2 = rbind(out[1,], rDNA)
rDNA[,2:ncol(rDNA)] = as.list(colSums(out2[,2:ncol(out2)]))
rDNA

out3 = rbind(out, rDNA)[c(1,2,4,3,5),]

tab_df(out3,
       col.header = names(out3),
       title = "Table 1. Summary statistics of segmental duplications in T2T CHM13 and GRCh38.",
       footnote = "peri: within 5 Mbp of the centromere; telo: within 500 kbp of the telomere; acro: within the short arms of the acrocentric chromosomes",
       show.footnote = TRUE, 
       show.rownames = FALSE )#sort.column = -3 )#, file="Tables/table1.html")

#webshot("Tables/table1.html", "Tables/table1.pdf", vwidth = 1200, vheight = 200)
#webshot("Tables/table1.html", "Tables/table1.png", vwidth = 1100, vheight = 200)
#write.xlsx2(x = data.table(out), file = "Tables/all_tables.xlsx", 
#            sheetName = "Table 1 SD Summary", row.names = FALSE)



achro_to_other = SEDEF[chr %in% ACHRO & ! chr2 %in% ACHRO, c("chr2","start2","end2")]
dim(achro_to_other)
sum(width(GenomicRanges::reduce(toGRanges(achro_to_other))))


#
#
# greater CN africans duplicons
#
#
rgn.df.with.genes = add_genes(rgn.df, all=TRUE)
rgn.df.with.genes$gene= rgn.df.with.genes$V4
sort((rgn.df.with.genes$gene))
unique(sort((rgn.df.with.genes$gene)))

genes = c("NBPF15", "GOLGA6L4", "NPIPB4","TBC1D3")

keep = grepl(paste(genes,collapse = "|"), rgn.df.with.genes$gene)
#clear_winner = c("chr9_3_12_cn_10", "chr1_5_2_cn_10", "chr13_1_79_cn_20");
#clear_winner = c(unique(rgn.df.with.genes$name[keep]),"chr6_161865995_161994506")
gene_rgn = unique(rgn.df.with.genes[keep,c("name","gene")], by="gene", fromLast=F);gene_rgn
gene_rgn = unique(rgn.df.with.genes[keep,c("name","gene","median")][order(median)], by="gene", fromLast=T) ; gene_rgn
#gene_rgn = rbind(gene_rgn, data.table(name="chr6_161865995_161994506", gene="LPA",median=30)); gene_rgn

winners = merge(rdg, gene_rgn , by="name")
#winners$genef = gsub("_.*", "",winners$gene)
#winners$genef=""
#lapply(genes,  function(x) winners$genef[grepl(x, winners$gene)] = x)
#unique(winners$genef)

make_gene_fam = c("NPIPB4"="NPIP", "TBC1D3K"="TBC1D3", "NBPF15"="NBPF",
                  "GOLGA6L4"="GOLGA", "SMA"="SMA")
winners$genef = make_gene_fam[winners$gene]
#winners$gene[str_starts(winners$gene,"GOLGA")]="GOLGA"
winners$label = paste(winners$genef, 
                      paste(paste(winners$chr, comma(winners$start), sep=":"), comma(winners$end), sep="-"),
                      sep="\n")
winners=winners[order(status)]
winners$label = factor(winners$label, levels = unique(winners$label))
winners$African = winners$super_pop=="AFR"

stats = winners %>% group_by(genef) %>%
  summarise(m = list(wilcox.test(CN[African], CN[!African], alternative="g")), 
            `AFR mean`= mean(CN[African]),
            `non-AFR mean`= mean(CN[!African])) %>%
  rowwise() %>%
  mutate(`p value` = m$p.value, W = m$statistic, method=m$method, alternative = m$alternative)

stats$`p value` = as.character(signif(stats$`p value`,digits = 3))

core_afr_cn.df = data.table(stats)[,-c("m","method","alternative", "W")]
tab_df(core_afr_cn.df, 
       title = "The copy number of core duplicons in Africans vs non-Africans",
       footnote = paste(unique(stats$method),",",unique(stats$alternative)), show.footnote = T,
       col.header = c("gene family", "AFR CN mean", "non-AFR CN mean", "p-value") )


#write.xlsx2(x = core_afr_cn.df, file = "Tables/all_tables.xlsx", sheetName = "Core duplicons in AFR vs non-AFR",
#           row.names=FALSE, append=TRUE)



#
#
#  Assembly info
#
#
#
stats.df = fread("../assemblies_for_anlysis/sample_info/sd_freeze_asm_stats.tbl")
asm.df = merge(fread("../assemblies_for_anlysis/sample_info/Master_SD_freeze.tbl"),stats.df, by.x="fasta", by.y="file") 
asm.df[sample =="GRCh38chrOnly"]$sample="GRCh38"



asm.df$Gender = as.character(asm.df$Gender)
asm.df[is.na(Gender)]$Gender = ""
asm.df[Gender == "2"]$Gender = "Female"
asm.df[Gender == "1"]$Gender = "Male"


asm.df[sample == "HG002", c("SuperPop", "Population", "Gender")] = data.table(SuperPop="EUR", Pop="AJ", Gender="Male")
asm.df[sample == "CHM13" | sample == "CHM1", c("SuperPop", "Population", "Gender")] = data.table(SuperPop="EUR", Pop="NA", Gender="Female")
asm.df[sample == "GRCh38", c("SuperPop", "Population", "Gender")] = data.table(SuperPop="NA", Pop="NA", Gender="NA")

asm.df$Species = "Human"
asm.df[grepl("_", asm.df$sample), c("sample", "Species")]  = asm.df[grepl("_",sample), c("sample")] %>%separate(sample, into=c("sample", "Species"), "_")

asm.df$Species = factor(asm.df$Species, levels = unique(c("Human", asm.df$Species)) )
asm.df[Species != "Human"]$Gender = "Female"
asm.df[sample == "Clint"]$Gender = "Male"
asm.df[Species != "Human", c("SuperPop", "Population")] = NA


asm.tbl = asm.df[order(sample!="CHM13", sample!="GRCh38", Species, sample),c("sample", "hap", "Species", "Gender", "SuperPop","Population", "totalBp", "N50", "auN", "nSeqs")]
asm.tbl


tab_df(asm.tbl,
       title = "Table 2. Assemblies used for analysis",
       footnote = "",
       show.footnote = TRUE, 
       show.rownames = TRUE,
       file="Tables/table_assemblies.html")


#webshot("Tables/table_assemblies.html", "Tables/table_assemblies.pdf")
#webshot("Tables/table_assemblies.html", "Tables/table_assemblies.png")
#write.xlsx2(x = asm.tbl, file = "Tables/all_tables.xlsx", sheetName = "Assemblies used in analysis",
#           row.names=FALSE, append=TRUE)




#
#
# Loci of interest 
#
#
pull.df = fread("../sd_regions_in_hifi_wga/pull_by_regions_snake_results/pull_sd_regions/results.tbl")
pull.df$N = nrow(unique(pull.df[Species=="Human",c("sm","hap")])) - 1 # for hg38
pull.df[Species != "Human"]$N = 2


fs = Sys.glob("../sd_regions_in_hifi_wga/pull_by_regions_snake_results/pull_sd_regions/Minigraph/*.fai")
rgns = gsub(".all.fasta.fai", "", basename(fs))
l <- lapply(fs, fread, sep="\t")
names(l)=rgns

het.df <- bind_rows(l, .id="Region") %>% 
  separate(V1, c("sm","hap","count"), sep = "\\.|__")%>% 
  mutate(count = as.numeric(count)) %>%
  separate(sm,c("sm", "Species"), sep = "_") %>%
  tidyr::replace_na(list(Species="Human")); het.df


pull.df = merge(pull.df, het.df, all.x=T)

pull.df = pull.df[ !(Region == "NOTCH2NL" & length < 50000 & sm == "HG03125") ]

pull.tbl = pull.df %>% filter(Region!="SRGAP2B_D" & Region != "SRGAP2C" & Region != "SRGAP2A" & Species=="Human") %>%
  group_by(Region) %>%
  summarise(`CHM13 (kbp)` = comma(length[sm=="CHM13"]/1000),
            `GRCh38 (kbp)` = comma(length[grepl("GRCh38", sm)]/1000),
            `Average length (kbp)` = comma(mean(length)/1000),
            `s.d. (kbp)` = comma(sd(length/1000)),
            `# resolved haplotypes`= n()-1, # for hg38
            `% resolved` = 100*(`# resolved haplotypes`)/unique(N),
            `% heterozygous haplotypes` = 100*sum(!is.na(V5))/(`# resolved haplotypes`+1)
            )

be_coords = fread("../sd_regions_in_hifi_wga/pull_by_regions_snake_results/regions.bed", col.names = c("chr", "start","end","Region"))

pull.tbl=merge(pull.tbl, be_coords, all.x=T)

tab_df(pull.tbl,
       title = "Table 2. Assemblies of evolutionary and biomedically important loci",
       col.header = names(pull.tbl),
       footnote = "Resolved haplotypes extend 50 kbp into unique sequence",
       show.footnote = TRUE, 
       show.rownames = TRUE
       )
#
#
# SVs in the evolutionary loci
#
#
sv = fread("../sd_regions_in_hifi_wga/pull_by_regions_snake_results/pull_sd_regions/sv.results.tbl")
colnames(sv)[1:6] = c("contig", "SV.start", "SV.end", "bubble.contig", "bubble.start", "bubble.end")
sv.tab = sv %>% filter(description != "Syntenic") %>%
  filter(!(region == "NOTCH2NL" & length < 50000 & contig == "HG03125.mat__1")) %>%
  separate(description, into = c("SV.type", "SV.kbp"), sep = "_|(kbp)") %>%
  mutate(SV.kbp = as.numeric(SV.kbp)) %>% filter(SV.kbp > 0) %>%
  separate(contig, into = c("contig", "# of haps"), sep='__')  %>%
  group_by(region, length, contig, `# of haps`, haplotypes, SV.type) %>%
  summarise(`#` = n(), `kbp` = paste(SV.kbp, collapse=";"), `total_kbp` = sum(SV.kbp)) %>%
  pivot_wider(names_from = `SV.type`, values_from = c("#", "kbp", "total_kbp"),
              names_sort=T) %>%
  mutate(`All SV kbp`= sum(total_kbp_DEL,total_kbp_INS, total_kbp_INV, na.rm=T)) %>%
  mutate(`Total # SVs`= sum(`#_DEL`,`#_INS`, `#_INV`, na.rm=T))%>%
  relocate(haplotypes, .after = last_col()) %>%
  relocate(`Total # SVs`, .after = `# of haps`)%>%
  relocate(`All SV kbp`, .after = `# of haps`) %>%
  replace_na(value = 0)
  

#
# dupmakser duplicons
# 
#
dm_diffs = dups[,1:7]


#
#
# genicn CN differences 
#
#
to_table = unique(rgn.df.with.genes[,c("gene","chr","start","end","ORF",
                                       "mean","median","std",
                                       "status","chm13_cn", "hg38_cn",
                                       "CHM13","GRCh38")], by=c("gene","chr","start","end"))[order(status,-ORF,gene)]

setnames(to_table, 
         old = c('CHM13','GRCh38', "chm13_cn", "hg38_cn"), 
         new = c('# with CHM13 CN','# with GRCh38 CN',
                 "CHM13 CN", "GRCh38 CN"), 
         skip_absent = T)
write.table(to_table, file="Tables/CN_of_non_syntenic_genes.tsv", sep="\t", row.names = F)


#
#
# # of SDs window counts (see sd_density_comparison.r )
#
#
sd_com = copy(sd_count[, c(1:3,6:10)]) %>% arrange(-diff)
names(sd_com) = c("chr", "start","end", '# CHM13 SDs', '# GRCh38 SDs', 'CHM13 SD bp', "GRCh38 SD bp", "# difference")

#
#
# peri to acro
#
#
peri_to_acro = SEDEF[peri & ( !chr %in% ACHRO) & acro2] %>% group_by(chr) %>%
  summarise(
            start = min(start),
            end = max(end),
            `# SDs` = n(),
            `Peri non-redudant kbp` = sum(width(GenomicRanges::reduce(toGRanges(chr,start,end))))/1e3, 
            `Redudant aligned kbp` = sum(matchB + mismatchB)/1e3, 
            `Ave. % ID` = mean(100*fracMatch)) %>%
  arrange(-`Redudant aligned kbp`); peri_to_acro

#
#
# rDNA collapses in the assembly 
#
#
readdepth = 32.67149444733489 
rdna_rd = fread("../sda_out/v1_collasped_regions_to_look_at/sd.collapsed.bed",
                col.names = c("chr","start","end","mean coverage","median coverage", "length"))
rdna_rd.tab = rdna_rd %>% filter(`mean coverage` > 100) %>%
  mutate(`Expected kbp of uncollapsed sequence` = `mean coverage` * length / (1e3 * readdepth) ) 

sum(rdna_rd.tab$`Expected kbp of uncollapsed sequence`)
 

#
#
# largest duplicons in the genomes
#
#


#Add the GC numbers
gc = DUPLICON_GC 

all_duplicons = rbind(DM, DM_38) %>% group_by(Assembly, duplicon, ancestral_position) %>%
  summarise(count=n(), kbp=sum(end-start)/1e3, `average length (kbp)` = mean(end-start)/1e3 ) %>%
  group_by(Assembly, duplicon, ancestral_position) %>%
  pivot_wider(names_from = Assembly, values_from = c(count, kbp, `average length (kbp)`), names_glue = "{Assembly} {.value}") %>%
  mutate(`Difference kbp` = `T2T CHM13 kbp` - `GRCh38 kbp`) %>%
  arrange(- `Difference kbp` ) %>%
  mutate(larger = case_when( `T2T CHM13 kbp` >= `GRCh38 kbp` ~ "T2T CHM13", 
                             T ~ "GRCh38")
         ) %>% 
  left_join(gc) %>%
  left_join(DM_GENES)

data.table(all_duplicons)
openxlsx::write.xlsx(
  list(`Largest duplicons`=all_duplicons),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Largest duplicons.xlsx")
#
#
# geneic duplicons with faimly expansions in the top 30 different
#
#
larger_duplicon_families = fread("../top_30_duplicons_gene_intersect/duplicon.gene.fam.expanded.tbl", 
                                 col.names = c("count", "Gene family")) %>%
  filter(count > 1)
  
  
openxlsx::write.xlsx(
  list(`Expanded gene families`=larger_duplicon_families),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Expanded_gene_families_within_the_30_most_expanded_duplicons.xlsx")





#
#
# rDNA reads
#
#
hifirdna = fread("../t2t_globus_share/team-rdna/reads-per-chromosome-20201120/rdna_read_stats.tbl") %>% 
  filter(type == "hifi") %>%
  mutate(`Expected bp of rDNA` = totalBp/readdepth )
hifirdna
sum(hifirdna$`Expected bp of rDNA`)/1e6

#
#
# write all output tables 
#
#


openxlsx::write.xlsx(
  list(`Table 1`=out3),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Table1.xlsx")

openxlsx::write.xlsx(
  list(`Core Duplicon CN in AFR vs non-AFR`=core_afr_cn.df),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Core_duplicon_CN_in_AFR_vs_non_AFR.xlsx") 

openxlsx::write.xlsx(
  list(`Evolutionary and biomedical loci`=pull.tbl),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Evolutionary_and_biomedical_loci.xlsx") 

openxlsx::write.xlsx(
  list(`Evolutionary and biomedical loci SVs`=sv.tab),
  headerStyle = hs,numFmt = "COMMA", colWidths=80, gridLines=TRUE, colNames = TRUE, 
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Evolutionary_and_biomedical_loci_SVs.xlsx") 


openxlsx::write.xlsx(
  list(`Assemblies used for analysis`=asm.tbl),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Assemblies_used_for_analysis.xlsx") 

openxlsx::write.xlsx(
  list( `New genes by liftoff`=data.table(NEW_GENES)[cdslen>=200]),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/New_genes_by_liftoff.xlsx") 

openxlsx::write.xlsx(
  list(`CN of non-syntenic genes`=to_table),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/CN_of_non_syntenic_genes.xlsx") 

openxlsx::write.xlsx(
  list(`Non-sytenic regions`= NEW),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Non_sytenic_regions.xlsx") 

openxlsx::write.xlsx(
  list(`CN in all genotyped regions` = rgn.df),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/CN_in_all_genotyped_regions.xlsx") 

openxlsx::write.xlsx(
  list(`Number of SDs per window` = sd_com),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Number_of_SDs_per_window.xlsx") 


openxlsx::write.xlsx(
  list(`Differences in duplicon content` = dm_diffs),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Differences_in_duplicon_content.xlsx") 

openxlsx::write.xlsx(
  list(`Genic SD expansions in T2T relative to Clint` = PTR_MIS),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Genic_SD_expansions_in_T2T_relative_to_Clint.xlsx") 


openxlsx::write.xlsx(
  list(`SDs from pericentromeric regions to acrocentric arms` = peri_to_acro),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Pericentromeric_regions_to_acrocentric_arms.xlsx") 

openxlsx::write.xlsx(
  list(`rDNA coverage collapses` = rdna_rd.tab),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/rDNA_coverage_collapses.xlsx") 

openxlsx::write.xlsx(
  list(`Chr specific rDNA read stats` = rdna_rd.tab),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/rDNA_coverage_collapses.xlsx") 


#
# large tables 
#
openxlsx::write.xlsx(
  list(`CN of all SGDP samples` = rdg),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Large Tables/CN_of_all_SGDP_samples.xlsx") 

openxlsx::write.xlsx(
  list(`Segmental duplications` = SEDEF),
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Large Tables/Segmental_duplications.xlsx") 


#
#
# divergence summary table  
#
#
sty2 = createStyle(numFmt="0.00")
openxlsx::write.xlsx(
  list(`Divergence`=divergence_summary),
  style = sty2,
  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
  "~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables/Divergence.xlsx")




