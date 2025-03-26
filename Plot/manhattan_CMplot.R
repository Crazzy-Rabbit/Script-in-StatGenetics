#! /public/home/shilulu/Lu_mamba/Configs/envs/R4.2.0/bin/Rscript

#  file columns should has SNP CHR BP P
args        <- commandArgs(TRUE)
ptype       <- args[1]
infile      <- args[2]
outname     <- args[3]

library(CMplot)
library(data.table)
library(dplyr)

data        <- fread(infile) %>% filter(complete.cases(.))
if (!all(c("SNP", "CHR", "BP", "P") %in% names(data))) stop("The file you provide not contains columns SNP, CHR, BP, and P")
pltdt       <- subset(data, select=c(SNP, CHR, BP, P))
plt_filter  <- na.omit(pltdt[pltdt$P > 0 & pltdt$P < 1, ])  # remove na and keep 0 < p < 1
my_color    <- c('#1B2C62', '#4695BC')

manhattan   <- function (plt, out, ylim=NULL, dpi=500, height=5, width=15, highlight=NULL, highlight.pch=19) {
  CMplot(plt, type="p", plot.type="m", LOG10=TRUE, ylim=ylim,
         file="jpg", file.name=out, file.output=TRUE, dpi=dpi, 
         threshold=5e-8, threshold.col="red", threshold.lwd=2, threshold.lty=2, 
         highlight=highlight, highlight.pch=highlight.pch, 
         col=my_color, band=0, 
         cex=0.6, signal.cex=0.6, height=height, width=width)
}

qqplot    <- function(plt, out, dpi=500, width=10, height=10){
  CMplot(plt, plot.type="q", LOG10=TRUE, file="jpg", file.name=outname, dpi=dpi, 
       file.name=out, verbose=TRUE, cex=0.8, signal.cex=0.8, width=width, height=height, 
       threshold.col="red", threshold.lty=2, conf.int=TRUE, conf.int.col=NULL, file.output=TRUE)
}

if (ptype == "manhattan") {
       if (any(-log10(plt_filter$P) > 50)){
              manhattan(plt=plt_filter, out=outname, ylim=c(0, 50))
       } else {
            manhattan(plt=plt_filter, out=outname)  
       }
} elseif (ptype == "qqplot") {
       qqplot(plt=plt_filter, out=outname)
}