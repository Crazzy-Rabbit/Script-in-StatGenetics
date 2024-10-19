# # # run_SMRinR_fiter.r
library(dplyr)
library(stringr)
library(data.table)

# r function run smr qtl to trait
run_smr_qtl <- function(gwas, qtl, outname){
    SMR="/public/home/lijunpeng/smr-1.3.1-linux-x86_64/smr-1.3.1"
    bfile="/public/home/lijunpeng/smr-1.3.1-linux-x86_64/g1000_eur/g1000_eur"

    system(paste(SMR, "--bfile", bfile, 
                      "--beqtl-summary", qtl, 
                      "--gwas-summary", gwas, 
                      "--diff-freq-prop 0.5", 
                      "--maf 0.01", 
                      "--cis-wind 2000", 
                      "--out", outname, 
                      "--thread-num 10"), intern=TRUE) -> mylog
    writeLines(mylog, paste0(outname, ".log"))
}

# p_SMR & p_HEIDI filter
run_filter < function(smr, outname){
    file = fread(smr)
    filter = file %>%
               filter(p_SMR <= 0.05/nrow(file) & p_HEIDI >= 0.01)

    write.table(filter_out, file=paste0(outname, "_", sprintf("%.2e", 0.05/nrow(file)), "_pass_smr_0.01_pass_heidi.txt"), 
            sep="\t", row.names=FALSE, quote=FALSE)
}


eQTL = "/public/home/gaikai/data/multi_tissue_xQTL/eqtl_all_tissue/Whole_Blood"
eQTLGen = "/public/home/gaikai/data/multi_tissue_xQTL/eqtl_all_tissue/Whole_Blood"
GWAS = "reformatMETAL.gz"

# select GWAS file cols: SNP A1 A2 freq beta SE N 
setwd("/public/home/shilulu/project_hearing-loss/new_run/all_meta/")
gwas = fread(GWAS)[, .(SNP, A1, A2, freq, beta, SE, N)]

# run smr and filter
setwd("~/SMR/whole_blood/")
outprx = str_replace(GWAS, ".gz", "_whole_blood")

for (qtl in c(eQTL, eQTLGen)){
    qtlname = ifelse(grepl("eQTL", qtl), "eQTL", "eQTLGen")
    outname = paste0(outprx, qtlname)
    run_smr_qtl(gwas, qtl, outname)
    run_filter(paste0(outname,"smr"), outname)
}
