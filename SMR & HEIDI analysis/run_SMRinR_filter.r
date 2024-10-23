#***********************************************#
# File   :   run_SMRinR_filter.r                 #
# Time   :   2024/10/23 20:50:51                #
# Author :   Lulu Shi                           #
# Mails  :   crazzy_rabbit@163.com              #
# link   :   https://github.com/Crazzy-Rabbit   #
#***********************************************#
library(dplyr)
library(optparse)
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
    smr_tsh = sprintf("%.2e", 0.05/nrow(smr_file))
    flt_file = smr_file %>% filter(p_SMR <= 0.05/n() & p_HEIDI >= 0.01)

    write.table(flt_file, file=paste0(outname, "_", smr_tsh, "smr_0.01heidi.txt"), sep="\t", row.names=FALSE, quote=FALSE)
}

# args <- commandArgs(TRUE)
# gwas <- args[which(args == "--gwas") + 1]
# qtl  <- args[which(args == "--qtl") + 1]
# out  <- args[which(args == "--out")  + 1]
# func <- ifelse("--func" %in% args, args[which(args == "--func") + 1], "all")

option_list <- list(
    make_option(c("--func", type="character"), default="both", help="Function to run: 'run_smr_qtl', 'run_filter', or 'both'"),
    make_option(c("--gwas"), type="character", help="GWAS summary data to run smr"),
    make_option(c("--qtl"), type="character", help="qtl file to run smr"),
    make_option(c("--out"), type="character", help="out prefix"),
    make_option(c("--smr"), type="character", help="SMR file for filtering")
)
opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

gwas <- opts$gwas
qtl  <- opts$qtl
out  <- opts$out
func <- opts$func
smr  <- opts$smr

if (func == "run_smr_qtl"){
    if (is.null(gwas) || is.null(qtl)) {
        stop("When func is 'run_smr_qtl', --gwas and --eqtl must be provided.")
    }
    run_smr_qtl(gwas, qtl, out)
} else if (func == "run_filter"){
    if (is.null(smr)) {
        stop("When func is 'run_filter', --smr must be provided")
    }
    run_filter(smr, out)
} else{
    if (is.null(gwas) || is.null(qtl)) {
        stop("When func is not provided, --gwas and --qtl must be provided.")
    }
    run_smr_qtl(gwas, qtl, out)
    run_filter(paste0(out, ".smr"), out)
}
