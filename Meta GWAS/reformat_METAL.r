library(dplyr)
library(data.table)

df <- fread("ARHL_EAS_MVP_BBJ1.tbl")

old_col = c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "TotalSampleSize")
new_col = c("SNP", "A1", "A2", "freq", "beta", "SE", "p", "N")
setnames(df, old_col, new_col)

df_filter <- df %>%
            filter(Indix==2) %>%
            mutate(CHR = Chr / 2) %>%
            mutate(BP = Pos /2) %>%
            select(CHR, BP, SNP, A1, A2, freq, beta, SE, p, N)

fwrite(df_filter, file="ARHL_EAS_MVP_BBJ_reformatMETAL.gz", sep="\t", row.names=FALSE)
