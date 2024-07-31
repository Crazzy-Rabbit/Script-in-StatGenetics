#! /public/home/shilulu/Lu_mamba/Configs/envs/R4.2.0/bin/R

## manhattan plot
#  file columns should has SNP CHR BP P
args <- commandArgs(TRUE)
infile <- args[1]
outname <- args[2]

my_color <- c('#1B2C62', '#4695BC')
data <- read.table(infile, header=TRUE)
pltdt <- subset(data, select=c(SNP, CHR, BP, P))
plt_filter <- pltdt[pltdt$P > 0, ]

library(CMplot)
CMplot(plt_filter, type="p", plot.type="m", LOG10=TRUE, file="jpg", file.name=outname, dpi=500,
       threshold=5e-8, threshold.col="red", threshold.lwd=2, threshold.lty=1, verbose=TRUE,
       col=my_color, band=0, cex=0.8, signal.cex=0.8, file.output=TRUE, height=5, width=15)
