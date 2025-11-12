#----------------------------------------------------#
#--------------        SuSiE      -------------------#
#----------------------------------------------------#
library(susieR);
library(Matrix);
library(data.table);

readEigenBin <- function(block_dir, info, block) {
    sinfo <- info[Block == block][order(Index)]
    if (nrow(sinfo) < 3) { cat(" skip: SNPs <", nrow(sinfo), "\n"); next }

    ld_file <- file.path(block_dir, sprintf("block%s.eigen.bin", block))
    con <- file(ld_file, "rb")
    m <- readBin(con, what=integer(), n=1, size=4)
    k <- readBin(con, what=integer(), n=1, size=4)
    sumPosEigVal <- readBin(con, what=numeric(), n=1, size=4)
    eigenCutoff  <- readBin(con, what=numeric(), n=1, size=4)
    evals <- readBin(con, what=numeric(), n=k, size=4)
    U <- matrix(readBin(con, what=numeric(), n=m*k, size=4), nrow=m, ncol=k)
    close(con)
    stopifnot(nrow(sinfo) == m)

    return(list(evals=evals, U=U, sinfo=sinfo))
}

lowrank_nearPD <- function(U, evals, sinfo, eps=1e-4) {
    R <- U %*% diag(evals, nrow=length(evals)) %*% t(U)
    R <- (R + t(R)) / 2
    diag_add <- pmax(eps, 1 - diag(R))
    R <- R + diag(diag_add)
    s <- 1 / sqrt(pmax(1e-10, diag(R)))
    R <- diag(s) %*% R %*% diag(s)
    R <- (R + t(R)) / 2
    if (min(eigen(R, symmetric=TRUE, only.values=TRUE)$values) < 1e-6) {
        R <- as.matrix(nearPD(R, corr=TRUE)$mat)
    }

    rownames(R) <- colnames(R) <- sinfo$ID
    return(R)
}

# ========= susieR with low-rank eigen-LD (master snp.info) =========
gwas      <- "All_MVP_Trpchevska_De-Angelis_BBJ.imputed.ma" 
ld_dir    <- "/public/home/shilulu/Wulab/GCTB/Eigen_decomposition_1M_ukbEUR_HM3"
snp_info  <- "snp.info"
out_dir   <- "susie"
dir.create(out_dir, showWarnings = FALSE)

#- read data
sumstats <- fread(gwas)
info     <- fread(snp_info)
idx      <- match(info$ID, sumstats$SNP)
if (anyNA(idx)) stop ("Some SNPs missing in summary data.")

sumstats <- sumstats[idx]
stopifnot(all(info$ID == sumstats$SNP))

#- flip alleles
flip <- info$A1 == sumstats$A2 & info$A2 == sumstats$A1
sumstats[flip, `:=`(freq=1-freq, b=-b)]

#- z scores
if (!"Z" %in% names(sumstats)) sumstats[, Z := b / se]

#- for each block
diag_list <- list()
cs_list   <- list()
pip_list  <- list()
for (blk in sort(unique(info$Block))) {
    cat("Block", blk, "...\n")

    eig <- readEigenBin(ld_dir, info, blk)
    evals <- eig$evals
    U     <- eig$U
    sinfo <- eig$sinfo

    #- match SNPs in gwas and eigen.bin
    ss_blk <- sumstats[SNP %in% sinfo$ID]
    ss_blk <- ss_blk[match(sinfo$ID, ss_blk$SNP)]
    if (nrow(ss_blk) < 3) { cat("  skip: matched SNPs <", nrow(ss_blk), "\n"); next }

    z <- ss_blk$Z
    N <- median(ss_blk$N, na.rm=TRUE)

    #- low-rank LD + residual diag + corrlation + nearPD
    R <- lowrank_nearPD(U, evals, sinfo)

    #- diagnostic of SuSiE model fitting
    lambda  <- estimate_s_rss(z=z, R=R, n=N)
    #- fine-mapping with SuSiE
    res <- susie_rss(z=z, R=R, n=N, L=10, estimate_residual_variance=TRUE, estimate_prior_variance=TRUE, check_R=FALSE)
    pip <- data.frame(SNP=names(res$pip), PIP=as.numeric(res$pip))
    cs_list <- susie_get_cs(res)
    cs <- data.frame(CS_ID = rep(names(cs$cs), lengths(cs$cs)),
                    SNP    = unlist(lapply(cs$cs, function(x) names(pip)[x])),
                    PIP    = unlist(lapply(cs$cs, function(x) pip[x])),
                    Coverage = rep(cs$coverage, lengths(cs$cs)))
    
    diag_list[[as.character(blk)]] <- lambda
    cs_list[[as.character(blk)]]   <- cs
    pip_list[[as.character(blk)]]  <- pip
}

#- save results
fwrite(rbindlist(diag_list, idcol="Block"), file=file.path(out_dir, "susie.lambda"), sep="\t")
fwrite(rbindlist(cs_list, idcol="Block"), file=file.path(out_dir, "susie.cs"), sep="\t")
fwrite(rbindlist(pip_list, idcol="Block"), file=file.path(out_dir, "susie.pip"), sep="\t")


#--------------------------------
# result plot SuSiE
#--------------------------------
#- plot diagnosis hist
setwd("F:/Github/PHD_job/2_project_hearing loss/NC_revision/gctb/gwfm/SuSiE")
library(data.table)
library(ggplot2)

dat = fread("arhl.lambda")
p = ggplot(dat, aes(x=V1)) +
  geom_histogram(bins=30, fill="grey70", color="black") +
  labs(x=expression(lambda), y="Count")+
  theme_bw()+
  theme(axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.line = element_line(color="black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(color=NA))
ggsave("susie_diagnostic_hist.png", p, width=8, height=4, dpi=600)


# ========= susieR with reference ld such 1KG =========
#- generate ld matrix with plink
plink --bfile 1KG --r square --extract snp.list --write-snplist --out gwas.snps

#- susieR
gwas      <- "All_MVP_Trpchevska_De-Angelis_BBJ.imputed.ma" 
ldprefix  <- "gwas.snps"
out_dir   <- "susie"
dir.create(out_dir, showWarnings = FALSE)

library(susieR)
library(data.table)
library(Matrix)

sumstats <- fread(gwas)
info     <- fread(paste0(ldprefix, ".snplist"))
R        <- fread(paste0(ldprefix, ".ld"))
R        <- as.matrix(R)

idx      <- match(info$ID, sumstats$SNP)
sumstats <- sumstats[idx]
stopifnot(all(info$ID == sumstats$SNP))

z        <- sumstats$b / sumstats$se
n        <- median(sumstats$N, na.rm=TRUE)
res      <- susie_rss(z, R, L=10, estimate_residual_variance=T, estimate_prior_variance=TRUE, check_R=FALSE, n=n)
pip      <- data.frame(SNP=names(res$pip), PIP=as.numeric(res$pip))
cs_list  <- susie_get_cs(res)
cs       <- data.frame(CS_ID     = rep(names(cs$cs), lengths(cs$cs)),
                        SNP      = unlist(lapply(cs$cs, function(x) names(pip)[x])),
                        PIP      = unlist(lapply(cs$cs, function(x) pip[x])),
                        Coverage = rep(cs$coverage, lengths(cs$cs)))


#- plot z score comparison
susie_plot(z, y = "z", b=b)

#- plot PIP comparison
susie_plot(res, y="PIP", b=b)