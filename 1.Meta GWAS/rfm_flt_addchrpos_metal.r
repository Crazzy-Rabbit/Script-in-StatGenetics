library(dplyr)
library(data.table)

# reformat the METAL out file
# filter the SNP which only occured at one study
METAL_rfm_flt <- function(file){
  df <- fread(paste0(file, "1.tbl.gz"))

  old_col = c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "TotalSampleSize")
  new_col = c("SNP", "A1", "A2", "freq", "beta", "SE", "p", "N")
  setnames(df, old_col, new_col)

  df_rfm <- df %>%
    filter(Indix > 1) %>%
    mutate(A1 = str_to_upper(A1), A2 = str_to_upper(A2))
  out_col <- df_rfm[, .(SNP, A1, A2, freq, beta, SE, p, N)]
  return(out_col)
}

# add chtr and pos information
add_chr_pos <- function(df) {
  setwd("/public/share/wchirdzhq2022/Wulab_share/dbSNP/GRCh37")
  db = fread(paste0(chr, ".txt"), header=FALSE, col.names=c("CHR", "POS", "SNP"))
  df_merge = merge(df, db, by="SNP", all=FALSE)
  return(df_merge)
}

setwd("/public/home/shilulu/project_hearing-loss/new_run")
metal_file = c("ARHL_MVP_AJHG_BBJ")

flt_metal = METAL_rfm_flt(metal_file)
setwd("/public/home/shilulu/project_hearing-loss/new_run/all_meta")
fwrite(flt_metal, file=paste0(metal_file, "_reformatMETAL"), sep="\t", row.names=FALSE, quote=FALSE)
out_added = data.table()
for (chr in c(1:22)){
  out_added = rbind(out_added, add_chr_pos(flt_metal))
}
fwrite(out_added, file=paste0(metal_file, "_reformatMETAL_addchr.gz"), sep="\t", row.names=FALSE, quote=FALSE)
