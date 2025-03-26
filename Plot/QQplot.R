#! /public/home/shilulu/Lu_mamba/Configs/envs/R4.2.0/bin/R

args <- commandArgs(TRUE)
infile <- args[1]
outname <- args[2]

data <- read.table(infile, header=TRUE)
pltdt <- subset(data, select=c(SNP, CHR, BP, P))
plt_filter <- pltdt[pltdt$P > 0, ]
plt_filter$P <- -log10(plt_filter$P)

library(CMplot)
CMplot(plt_filter, plot.type="q", LOG10=FALSE, file="jpg", file.name=outname, dpi=500, 
       verbose=TRUE, cex=0.8, signal.cex=0.8, width=10, height=10, 
       threshold.col="red", threshold.lty=2, conf.int=TRUE, conf.int.col=NULL, file.output=TRUE)
