##########################
###     GWAS treat     ###
##########################
library(data.table);library(dplyr);library(tidyr)
#--------------------------------------------#
# QC for GWAS sumary statistic               #
#--------------------------------------------#
sums_col_qc <- function(gwas, SNP, A1, A2, freq, P,
                            CHR = NULL, POS = NULL,
                            BETA = NULL, SE = NULL,
                            OR = NULL, Z = NULL, N=NULL, MAF=NULL) {
    gwas_col = c(
        CHR   =   CHR,
        POS   =   POS, 
        SNP   =   SNP,
        A1    =    A1, 
        A2    =    A2,
        freq  =  freq,
        BETA  =  BETA,
        SE    =    SE,
        OR    =    OR,
        Z     =     Z,
        N     =     N,
        P     =     P
    )
    gwas_col = gwas_col[!sapply(gwas_col, is.null)]

    missing_required = setdiff(c("SNP", "A1", "A2", "freq", "P"), names(gwas_col))
    if (length(missing_required) > 0) {
        stop("Missing required columns: ", paste(missing_required, collapse = ", "))
    }
    has_b_se = all(c("BETA", "SE") %in% names(gwas_col))
    has_or   = "OR" %in% names(gwas_col)
    has_z    = "Z" %in% names(gwas_col)
    if (!(has_b_se || has_or || has_z)) {
        stop("Must provide at least one of the following: (BETA + SE), OR, or Z.")
    }

    gwas = fread(gwas, select = unname(gwas_col))
    setnames(gwas, old = unname(gwas_col), new = names(gwas_col))
    
    valid_p <- !is.na(gwas$P) & gwas$P >= 0 & gwas$P <= 1
    if (has_b_se){
        gwas = gwas[valid_p & !is.na(BETA) & !is.na(SE) & BETA != 0, ]
    } else if (has_or){
        gwas = gwas[valid_p & !is.na(OR), ]
    } else if (has_z){
        gwas = gwas[valid_p & !is.na(Z), ]
    }
    
    if (!is.null(MAF)){
        gwas = gwas[freq > MAF & freq < 1-MAF, ]
    }

    return(gwas)
}

#--------------------------------------------#
# N effect and reascale beta se              #
#--------------------------------------------#
rescale_b_se <- function(gwas, case, control){
  names(gwas) = c("SNP","A1","A2","freq","b","se","P","N")
  # n           = gwas$N
  # i           = dnorm(qnorm(1-K))/K
  # n_eff       = i^2*v*(1-v)/(1-K)^2 * n # Yang et al. Gen Epi 2010
  n_eff       = as.integer(4*case*control/(case + control))
  p           = gwas$freq
  z           = gwas$b / gwas$se
  se          = 1/sqrt(2*p*(1-p)*(n_eff+z^2))
  b           = z*se
  gwas$b      = b
  gwas$se     = se
  gwas$N      = n_eff
  
  return(gwas)
}

##########################
###  ploting function  ###
##########################
library(ggplot2);library(data.table);
library(dplyr);library(RColorBrewer)
#--------------------------------------------#
# plot manhattan and qq for GWAS summary     #
#--------------------------------------------#
# ex: png("test.png", width=2400, height=2400, res=500, type="cairo")
# qqplot(pvals)
# dev.off()
qqplot = function(pval, ylim=0, lab=1.4, axis=1.2){
  par(mgp=c(5,1,0))
  p1 <- pval 
  p2 <- sort(p1)
  n  <- length(p2)
  k  <- c(1:n)
  alpha <- 0.05
  lower <- qbeta(alpha/2, k, n+1-k)
  upper <- qbeta((1-alpha/2), k, n+1-k)
  expect <- (k-0.05)/n 
  biggest <- ceiling(max(-log10(p2), -log10(expect)))
  xlim <- max (-log10(expect)+0.1);
  if (ylim==0) ylim=biggest;
  plot(-log10(expect), -log10(p2), xlim=c(0, xlim), ylim=c(0, ylim),
       ylab=expression(paste("Observed ", "-", log[10], "(", italic(P), "-value)", sep="")),
       xlab=expression(paste("Expected ", "-", log[10], "(", italic(P), "-value)", sep="")),
       type="n", mgp=c(2,0.5,0), tcl=-0.3, bty="n", cex.lab=lab, cex.axis=axis)
  polygon(c(-log10(expect), rev(-log10(expect))), c(-log10(upper), rev(-log10(lower))),
          col=adjustcolor("grey", alpha.f=0.5), border=NA)
  abline(0,1,col="white", lwd=2)
  points(-log10(expect), -log10(p2), pch=20, cex=0.6, col=2)
}


#--------------------------------------------#
# plot manhattan for SMR results             #   
#--------------------------------------------#
# ex: eQTL = plot_SMRmhn("eqtl.smr", showtext=TRUE)
plot_SMRmhn <- function(smr, showtext=FALSE){
    dt = fread(smr)[, .(Gene, ProbeChr, Probe_bp, p_SMR, b_SMR, p_HEIDI)]
    plt = dt[,.(Gene, ProbeChr, Probe_bp, p_SMR)]
    colnames(plt) = c("Gene", "CHR", "BP", "P")

    chr_lengths <- plt %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP))
    data <- plt %>%
    arrange(CHR, BP) %>%
    mutate(chr_cumsum = cumsum(c(0, head(chr_lengths$chr_len, -1)))[CHR]) %>%
    mutate(BPcum = BP + chr_cumsum)
    axis_set <- data %>%
    group_by(CHR) %>%
    summarise(center = (max(BPcum) + min(BPcum)) / 2)
    
    threshold = 0.05 / nrow(dt)
    Gene = dt[p_SMR < threshold & p_HEIDI > 0.01, ]
    Gene$color = Gene[, ifelse(b_SMR > 0, brewer.pal(8,"Dark2")[1:2][2], brewer.pal(8,"Dark2")[1:2][1])]
    Gene$pch = Gene[, ifelse(b_SMR > 0, 24, 25)]
    hlt_dt = merge(Gene, data, by.x=c("Gene", "Probe_bp"), by.y=c("Gene", "BP"))
    hlt_dt$pch = as.factor(hlt_dt$pch)

    p = ggplot() +
    geom_point(data = data, aes(x = BPcum, y = -log10(P), color = factor(CHR)), alpha = 0.9) +
    geom_point(data=hlt_dt, aes(x=BPcum, y=-log10(P), shape=pch),fill=hlt_dt$color, color=hlt_dt$color, alpha=0.9, size=2.5) +
    geom_hline(yintercept=-log10(threshold), color="grey68", linetype="dashed") +
    scale_color_manual(values=rep(c("#bfbfbf", "#b3dfe7"), 22)) + 
    scale_shape_manual(values = c(24, 25)) + 
    scale_x_continuous(label=axis_set$CHR, breaks=axis_set$center, expand=expansion(mult = c(0, 0.01))) +
    scale_y_continuous(expand=expansion(mult = c(0, 0.01))) + 
    labs(x="Chromosome", y=expression(-log[10] (italic(P_SMR)))) +
    theme_classic() +
    theme(legend.position="none",
        axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=15, margin=margin(t=10), color="black"),
        axis.title.y = element_text(size=15, margin=margin(t=10), color="black"),
        axis.line = element_line(color="black"), 
        axis.ticks = element_line(color="black"))
    
    if (showtext) {
        p = p + geom_text_repel(
        data = hlt_dt,
        aes(x=BPcum, y=-log10(P), label = Gene),
        color=hlt_dt$color,
        size = 3,
        fontface="italic",
        max.overlaps = Inf,
        max.iter = 10000, 
        force = 2,
        force_pull = 0.01,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey50") 
    }
    return(p)
}
