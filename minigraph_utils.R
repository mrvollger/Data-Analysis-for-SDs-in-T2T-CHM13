CHRS <<- c(paste0("chr",seq(1,22)),"chrX", "chrY", "chrM", "chrMT")
NOYM = CHRS[! CHRS %in% c("chrY","chrMT","chrM")]
NOM = CHRS[! CHRS %in% c("chrMT","chrM")]
#
# colors and gene name
#
GRAY = "#2F4F4F"	
RED = "#af0404"
BLUE = "#3282b8"
NEWCOLOR = RED
OLDCOLOR = GRAY 



tri_bed <- function(f, s=.2, allowed_names=c()){
  #df = fread(glue("../sd_regions_in_hifi_wga/lpa/minimiro/temp_minimiro/{gene}_query.fasta.duplicons.extra")); df
  df=fread(f)
  names(df)[1:3]=c("chr","start","end")
  if(length(allowed_names) > 0){
    df = df[chr %in% allowed_names]
  }
  df = df[order(chr,start)]
  df$chr = factor(df$chr)
  df$y = as.numeric(df$chr)
  df$y = 1
  df$tri_id = 1:nrow(df)
  if("orient" %in% names(df)){
    df$strand = "+"
    df$strand[df$orient == "R"] = "-"
  } else {
    df$strand = df$V6
    df$color = sapply(strsplit(df$V9, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
  }
  # allow to be faceted together
  df$q = df$chr
  zs = s/3
  data.table(df %>% 
               rowwise() %>% 
               mutate(xs=list(c(start, start, end, end)), 
                      ys = case_when(
                        strand == "+" ~ list(c(y+s,y-s,y-zs,y+zs)),
                        strand == "-" ~ list(c(y-zs,y+zs,y+s,y-s))
                      )) %>% unnest(cols=c("xs","ys")))
}

# make the seq names more human readable

clean_names <- function(list){
  list = gsub(" X 1$", "", 
              gsub("__", " X ", list)
  )
  list = gsub("CHM13.pri", "CHM13", list)
  list = gsub("GRCh38chrOnly.pri", "GRCh38", list)
  return(list)
}

clean_df <- function(df, allowed){
  allowed = factor(unique(c("CHM13", "GRCh38",clean_names(allowed))),
                   levels = unique(c("CHM13", "GRCh38",clean_names(allowed))))
  for( x in c("r", "q", "chr", "Label")){
    if( x %in% names(df) ){
      # clean the names
      df[[x]] = clean_names(df[[x]])
      # remove extra names
      df = df[ df[[x]] %in% allowed]
      # set the factors
      df[[x]] = factor(df[[x]], levels = allowed)
    }
  }
  return(df)
}
