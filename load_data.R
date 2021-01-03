#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plotutils.R")

if(F){
  NS = readbed("../Assembly_analysis/Masked/chm13.draft_v1.0_plus38Y.Ns.bed", "Ns")
  FAI = FAI[order(chr)]
  CENS=readbed(glue("../Assembly_analysis/SEDEF/{V}.cen.bed"), "Centromere")
  CENS$gieStain = "gpos100"#"gvar" #acen" , "gneg"
  CYTO.df = data.table(FAI$chr, 0, FAI$chrlen,gieStain="geng")
  CYTO2.df=CYTO.df
  CYTO2.df$gieStain[CYTO.df$V1 %in% ACHRO ] = "gpos75"#"stalk"
  GENOME=toGRanges(CYTO.df) 
  CYTO = c(GENOME,toGRanges(CENS))
  CYTO2 = c(toGRanges(CYTO2.df),toGRanges(CENS))
  
  #
  # reading in data data frames
  #
  synt=readbed(glue("../Assembly_analysis/snyteny/{V}.snyteny_1Mbp.bed"), "synt")
  NEW=grtodf(GenomicRanges::setdiff(GenomicRanges::reduce(GENOME),toGRanges(synt)))
  NEW=NEW[!NEW$chr %in% c("chrY", "chrMT","chrM", NA)]
  
  GENES=readbed(glue("../Assembly_analysis/Liftoff/{V}.orf_only.bed"), "GENES")
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
  SEDEF_37 = readbed("../Assembly_analysis/SEDEF/hg19.no_alt.SDs.bed", "GRCh37")
  SEDEF_CELARA = readbed("../Assembly_analysis/SEDEF/Celera_WGSA.SDs.bed", "Celera WGSA")
  
  NEW_GENES = read_excel(glue("../Assembly_analysis/Liftoff/{V}.liftoff.summary.xlsx"), sheet="NewCopiesLongestCDS")

  METH_SD_GENES = fread("../t2t_globus_share/team-epigenetics/20200727_methylation/v1.0_methylation/SD_analysis/sd.transcripts.10kb_methAG.bed.gz")
  
  
  save.image("~/Desktop/Rdata/plotutils.data")
} else {
  load("~/Desktop/Rdata/plotutils.data")
}


