library(data.table)
library(optparse)
option_list <- list(
    make_option(c("--gwas", type="character"), help="GWAS summary data"),
    make_option(c("--SNP"), type="character", help="SNP colnames"),
    make_option(c("--A1"), type="character", help="A1 colnames for effect allele"),
    make_option(c("--A2"), type="character", help="A1 colnames for another allele"),
    make_option(c("--freq"), type="character", help="frequency colnames for effect allele"),
    make_option(c("--beta"), type="character", help="beta colnames for effect allele"),
    make_option(c("--se"), type="character", help="se colnames"),
    make_option(c("--pval"), type="character", help="p colnames"),
    make_option(c("--num"), type="character", help="sample size colnames or sample size number"),
    make_option(c("--outprefix"), type="character", help="outprefix for your file")
)
opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

gwas <- opts$gwas
SNP  <- opts$SNP
A1   <- opts$A1
A2   <- opts$A2
freq <- opts$freq
beta <- opts$beta
se   <- opts$se
pval <- opts$pval
num  <- opts$num

# set col names
gwas_dt <- fread(gwas)
if (!is.na(as.numeric(num)) && !grepl("[a-zA-Z]", num)){
    gwas_dt$N <- as.numeric(num)
    old_col <- c(SNP, A1, A2, freq, beta, se, pval, "N")
    new_col <- c("SNP", "A1", "A2", "freq", "beta", "se", "p", "N")
    setnames(gwas_dt, old_col, new_col)
} else {
    if (num %in% colnames(gwas_dt)) {
        old_col <- c(SNP, A1, A2, freq, beta, se, pval, num)
        new_col <- c("SNP", "A1", "A2", "freq", "beta", "se", "p", "N")
        setnames(gwas_dt, old_col, new_col)
    } else {
        stop(paste("Column", num, "neither numeric nor found in the data"))
    }
}

fwrite(gwas_dt, file=outprefix, sep="\t")
