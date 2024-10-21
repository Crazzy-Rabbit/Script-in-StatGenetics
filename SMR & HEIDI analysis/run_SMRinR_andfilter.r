library(dplyr)
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
run_filter <- function(smr, outname){
    smr_file = fread(smr)
    flt_file = smr_file %>% filter(p_SMR <= 0.05/n() & p_HEIDI >= 0.01)

    write.table(flt_file, file=paste0(outname, "_", sprintf("%.2e", 0.05/nrow(smr_file)), "smr_0.01heidi.txt"), sep="\t", row.names=FALSE, quote=FALSE)
}

args <- commandArgs(TRUE)
gwas <- args[which(args == "--gwas") + 1]
qtl  <- args[which(args == "--qtl") + 1]
out  <- args[which(args == "--out")  + 1]

run_smr_qtl(gwas, qtl, out)
runfilter(paste0(out, ".smr"), out)
