setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")

if(F){
  
  NS = readbed("../Assembly_analysis/Masked/chm13.draft_v1.0_plus38Y.Ns.bed", "Ns")
  NS$gieStain = "stalk"
  FAI = FAI[order(chr)]
  CENS=readbed(glue("../Assembly_analysis/SEDEF/{V}.cen.bed"), "Centromere")
  CENS$gieStain = "gpos100"#"gvar" #acen" , "gneg"
  CYTO.df = data.table(FAI$chr, 0, FAI$chrlen,gieStain="geng")
  CYTO2.df=CYTO.df
  CYTO2.df$gieStain[CYTO.df$V1 %in% ACHRO ] = "gpos75"#"stalk"
  GENOME=toGRanges(CYTO.df) 
  CYTO = c(GENOME,toGRanges(CENS), toGRanges(NS))
  CYTO2 = c(toGRanges(CYTO2.df),toGRanges(CENS), toGRanges(NS))
  
  CYTO_COLORS = getCytobandColors()
  
  #
  # reading in data data frames
  #
  synt=readbed(glue("../Assembly_analysis/snyteny/{V}.snyteny_1Mbp.bed"), "synt")
  NEW=grtodf(GenomicRanges::setdiff(GenomicRanges::reduce(GENOME),toGRanges(synt)))
  NEW=NEW[!NEW$chr %in% c("chrY", "chrMT","chrM", NA)]
  
  GENES=readbed(glue("../Assembly_analysis/Liftoff/{V}.orf_only.bed"), "GENES")
  GENES$gene=GENES$V4
  RM = readbed("../Assembly_analysis/Masked/{V}_repeatmasker.out.bed", "T2T CHM13", rm=T)
  SAT = RM[RM$type == "Satellite"]
  
  
  DM_BED = readbed("../Assembly_analysis/Masked/{V}_dupmasker_colors.bed", "T2T CHM13")
  DM = readbed("../Assembly_analysis/Masked/{V}.duplicons.bed", "T2T CHM13")
  DM_38 = readbed("../Assembly_analysis/Masked/hg38.no_alt.duplicons.bed", "GRCh38")
  DM_37 = readbed("../Assembly_analysis/Masked/hg19.no_alt.duplicons.bed", "GRCh37")
  
  SEDEF = readbed("../Assembly_analysis/SEDEF/{V}.SDs.bed", "T2T CHM13")
  SEDEF = rgntag(SEDEF, DM, "Duplicon")
  LOW = readbed("../Assembly_analysis/SEDEF/{V}.SDs.lowid.bed", "T2T CHM13")
  
  ENRICHED = readbed("../Assembly_analysis/SEDEF/{V}.sedef.enriched.bed", "highsd")[,1:4]
  
  
  SEDEF_38 = readbed("../Assembly_analysis/SEDEF/hg38.chr_only.SDs.bed", "GRCh38", chrfilt=TRUE)
  SEDEF_38 = rgntag(SEDEF_38, DM_38, "Duplicon")
  
  SEDEF_37 = readbed("../Assembly_analysis/SEDEF/hg19.no_alt.SDs.bed", "GRCh37")
  SEDEF_37 = rgntag(SEDEF_37, DM_37, "Duplicon")
  
  SEDEF_CELARA = readbed("../Assembly_analysis/SEDEF/Celera_WGSA.SDs.bed", "Celera WGSA")
  
  NEW_GENES = read_excel(glue("../Assembly_analysis/Liftoff/{V}.liftoff.summary.xlsx"), sheet="NewCopiesLongestCDS")

  DM_GENES = fread("tmp.dup.gene", col.names = c("duplicon","gene"));
  
  
  ALL_ALN = readbed(glue("../Assembly_analysis/Align/{V}.split.bed"), "All 5kbp windows")
  #all = all[!grepl('chrY', all$chr)] # filter out chrY since we have no chrY, we now have a y so comment
  SDS_ALN = readbed(glue("../Assembly_analysis/Align/{V}.split.sd.bed"), "SDs")
  FLANK_ALN = readbed(glue("../Assembly_analysis/Align/{V}.split.sdflank.bed"), "SD flanks")
  NONSD_ALN = readbed(glue("../Assembly_analysis/Align/{V}.split.nosd.bed"), "Non SD")
  
  
  save.image("~/Desktop/Rdata/plotutils.data")
} else {
  load("~/Desktop/Rdata/plotutils.data")
}




