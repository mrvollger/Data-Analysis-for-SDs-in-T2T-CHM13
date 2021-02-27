setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(data.table)
library(glue)

# chromosome names
CHRS <<- c(paste0("chr",seq(1,22)),"chrX", "chrY", "chrM", "chrMT")
NOYM = CHRS[! CHRS %in% c("chrY","chrMT","chrM")]
NOM = CHRS[! CHRS %in% c("chrMT","chrM")]

# colors to use 
GRAY = "#2F4F4F"	
RED = "#af0404"
BLUE = "#3282b8"
NEWCOLOR = RED
OLDCOLOR = GRAY 
COLORS <<- c(`T2T CHM13`=NEWCOLOR, GRCh37=BLUE, GRCh38=GRAY, `Celera WGSA`="#ede682", `HG00733 pat`="#96bb7c",`HG00733 mat`="#ade498", `WGAC`="#000000")
V="chm13.draft_v1.0_plus38Y"
ACHRO <<- paste0("chr",c(13,14,15,21,22))
FAI <<- fread(glue("data/{V}.fasta.fai"),col.names = c("chr","chrlen","x","y","z"))
FAI$chr = factor(FAI$chr, levels =  c(CHRS, unique(FAI$chr[which(!FAI$chr %in% CHRS)]) ) , ordered = TRUE)

# output directories
FIGURE_DIR="~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Figures"
SUPP="~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Figures/Misc"
TABLES="~/Google Drive/My Drive/Vollger CHM13 T2T SDs 2020/Tables"
LOCAL_DATA="~/Desktop/Rdata"

#
# get the helper functions
#
source("plotutils.R")

#
# load large data
#
if(F){
  FAI_38 = fread("data/hg38.chr_only.fa.fai", col.names = c("chr","length", "x","y","z"))
  
  #
  NS = readbed("data/chm13.draft_v1.0_plus38Y.Ns.bed", "Ns")
  NS$gieStain = "stalk"
  FAI = FAI[order(chr)]
  CENS=readbed(glue("data/{V}.cen.bed"), "Centromere")
  CENS$gieStain = "gpos100"#"gvar" #acen" , "gneg"
  CYTO.df = data.table(FAI$chr, 0, FAI$chrlen,gieStain="geng")
  CYTO2.df=CYTO.df
  CYTO2.df$gieStain[CYTO.df$V1 %in% ACHRO ] = "gpos75"#"stalk"
  GENOME=toGRanges(CYTO.df) 
  CYTO = c(GENOME,toGRanges(CENS), toGRanges(NS))
  CYTO2 = c(toGRanges(CYTO2.df),toGRanges(CENS), toGRanges(NS))
  
  
  CYTO_COLORS = getCytobandColors()
  CYTOtmp = toGRanges(fread("data/chm13.CytoBandIdeo.v2.txt", col.names = c("chr","start","end", "name","gieStain")))
  CYTO38 = getCytobands(genome="hg38")
  CYTOFULL = c(CYTOtmp,  CYTO38[seqnames(CYTO38) == "chrY"])
  #plotKaryotype(genome=GENOME, cytobands = CYTOFULL)
  #show_col(getCytobandColors())
  
  #
  # reading in data data frames
  #
  synt=readbed(glue("data/{V}.snyteny_1Mbp.bed"), "synt")
  NEW=grtodf(GenomicRanges::setdiff(GenomicRanges::reduce(GENOME),toGRanges(synt)))
  NEW=NEW[!NEW$chr %in% c("chrY", "chrMT","chrM", NA)]
  
  GENES=readbed(glue("data/{V}.orf_only.bed"), "GENES")
  GENES$gene=GENES$V4
  ALL_GENES=readbed(glue("data/{V}.all.bed"), "ALL_GENES")
  ALL_GENES$gene=ALL_GENES$V4
  
  RM = readbed("data/{V}_repeatmasker.out.bed", "T2T CHM13", rm=T)
  SAT = RM[RM$type == "Satellite"]
  
  
  DM_BED = readbed("data/{V}_dupmasker_colors.bed", "T2T CHM13")
  DM = readbed("data/{V}.duplicons.bed", "T2T CHM13")
  DM_38 = readbed("data/hg38.chr_only.duplicons.bed", "GRCh38")

  SEDEF = readbed("data/{V}.SDs.bed", "T2T CHM13")
  NONR_SD = GenomicRanges::reduce(toGRanges(SEDEF))
  SEDEF = rgntag(SEDEF, DM, "Duplicon")
  LOW = readbed("data/{V}.SDs.lowid.bed", "T2T CHM13")
  
  ENRICHED = readbed("data/{V}.sedef.enriched.bed", "highsd")[,1:4]
  
  
  SEDEF_38 = readbed("data/hg38.chr_only.SDs.bed", "GRCh38", chrfilt=TRUE)
  SEDEF_38 = rgntag(SEDEF_38, DM_38, "Duplicon")
  SEDEF_CELARA = readbed("data/Celera_WGSA.SDs.bed", "Celera WGSA")
  
  NEW_GENES = read_excel(glue("data/{V}.liftoff.summary.xlsx"), sheet="NewCopiesLongestCDS")

  DM_GENES = fread("data/gene_duplicon_map.tbl", col.names = c("duplicon","gene"));
  
  
  ALL_ALN = readbed(glue("data/{V}.split.bed"), "All 5kbp windows")
  SDS_ALN = readbed(glue("data/{V}.split.sd.bed"), "SDs")
  FLANK_ALN = readbed(glue("data/{V}.split.sdflank.bed"), "SD flanks")
  NONSD_ALN = readbed(glue("data/{V}.split.nosd.bed"), "Non SD")
  

  #
  # RDG 
  #
  RDG = NA
  regions = "TRIMMED_WINDOWS"
  GMM=""
  for(pop in c("archaics", "hgdp", "kmer_ref", "t2t_dp", "t2t_dp_nhp", "nhp")){
    t2 = fread(glue("data/rdg/{regions}.{pop}.wssd{GMM}.genotypes.tab"))
    RDG  = merge(RDG, t2)
  }
  POPINFO=fread("data/rdg/hgdp_manifest.txt")
  
  #
  # probe mappings
  #
  PROBE_MAPPINGS = readbed("data/probe.mappings.bed", "probes")
  PROBE_RESULTS = fread("data/probe.results.csv")
  
  save.image(glue("{LOCAL_DATA}/plotutils.data"), compress = FALSE)
} else {
  load(glue("{LOCAL_DATA}/plotutils.data"))
}



  

if(F){
  METH_SD_GENES = fread(glue("data/sd.transcripts.10kb_methAG.bed.gz"))
  METH_SD_GENES_002 = fread(glue("data/sd.HG002transcripts.10kb_methAG.bed.gz"))
  METH_CLUSTERS = fread("data/t2t_chm13v1.0_SD_clustered_methylation.bed")
  METH_CLUSTERS_002 = fread("data/t2t_HG002_SD_clustered_methylation.bed")
  
  METH_JUST_GENES = readbed("data/sd.transcripts.and.meth.bed","z")
  CPG_DATA =  fread("data/scaled_binned_slop_genes.cpg.bed") 
  
  save(METH_SD_GENES, METH_CLUSTERS, METH_SD_GENES_002, METH_CLUSTERS_002, METH_JUST_GENES, CPG_DATA,
       file = glue("{LOCAL_DATA}/meth.data"))
}else if (F){
  load(glue("{LOCAL_DATA}/meth.data"))
}






