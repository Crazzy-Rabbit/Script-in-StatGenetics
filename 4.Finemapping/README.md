#### 1、GWFM

```
#! /bin/bash

gctb="~/software/gctb"
ldref="~/GCTB_ref/ukbEUR_Imputed"
gwas="~/GWAS.txt"
out_block="~/gwfm/impu_block/GWAS_impu"
out_impu="GWAS.imputed"
ldfile="~/GCTB_ref/ukbEUR_Imputed/LD_rsq05.ld.txt"
annofile="~/GCTB_ref/annot_baseline2.2.txt"
gwas_impu="GWAS.imputed.ma"
```
step1: impute summary data
```
for i in {1..591}; do
  $gctb --ldm-eigen ${ldref} --gwas-summary ${gwas} --impute-summary --block $i --thread 10 --out ${out_block}
done
```

step2: combined block file to one 
```
$gctb --gwas-summary GWAS_impu --merge-block-gwas-summary --out ${out_impu} 
```

step3: Running GWFM
```
outprx=$(basename -- ${gwas_impu} ".imputed.ma")

$gctb --sbayes RC --ldm-eigen $ldref --annot $annofile --gwas-summary $gwas_impu --n-dist-aut --write-mcmc-bin --thread 10 --out ${outprx}_SbayesRC
```

step4: credible sets
```
$gctb --cs --ld-file ${ldfile} --pip 0.9 --mcmc-samples $${outprx}_SbayesRC --out ${outprx}_SbayesRC_finemapping
```


#### 可视化结果
###### 单个变异
```
source("~/script/plot_GWFM.r")
gwas="~/GWAS_addchr.gz"
gwfm="~/gwfm/GWAS_SbayesRC.snpRes"
genelist="~/script/plot_smr/glist_hg19_refseq.txt"

snp = "rs1126809"

PData = ReadPvalueFromFiles(gwas=gwas, gwfm=gwfm, glist=genelist, windowsize=300000, highlight=snp)
pdf("rs1126809_test.pdf",width = 8,height = 8)
MultiPvalueLocusPlot(data=PData)
dev.off()
```
###### 多个变异循环
```
source("~/script/plot_GWFM.r")
gwas="~/GWAS_addchr.gz"
gwfm="~/gwfm/GWAS_SbayesRC.snpRes"
genelist="~/script/plot_smr/glist_hg19_refseq.txt"

SNP = fread("~/GWAS_SbayesRC_finemapping.gcs")
causal = SNP[PIP > 0.9, ][ ,SNP]

causal = data.table(
    rsID = c("rs62623452", "rs13147559", "rs2242416", "rs9493627", "rs35094336", "rs143282422", "rs1126809",
            "rs35887622", "rs2289702", "rs1805005", "rs118174674", "rs61734651", "rs45598239", "rs36062310"),
    CHR = c("2", "4", "6", "6", "8", "10", "11", "13", "15", "16", "18", "20", "21", "22")
)

dir_path = "~/gwfm/plot"
dir.create(dir_path, showWarnings = FALSE)

for (i in seq_len(nrow(causal))) {
    snp = causal[i]$rsID
    chr = causal[i]$CHR
    print(paste("processing: ", snp))
    outprex = file.path(dir_path, sprintf("chr%s_%s_gwfm.pdf", chr, snp))
    PData = ReadPvalueFromFiles(gwas=gwas, gwfm=gwfm, glist=genelist, windowsize=300000, highlight=snp)
    pdf(outprex,width = 8,height = 8)
    MultiPvalueLocusPlot(data=PData)
    dev.off()
}
```
###### 结果示例
![1762479992443](image/GCTB_GWFM/1762479992443.png)

#### 2、SuSiE

step 1: generate ld matrix with plink
```
plink --bfile 1KG --r square --extract snp.list --write-snplist --out gwas.snps
```

step 2: Run SuSiE
```
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

```
step 3: plot gwas locus
```
#- plot z score comparison
susie_plot(z, y = "z", b=b)
```
![1762916414358](image/README/1762916414358.png)

step 4: plot PIP
```
#- plot PIP comparison
susie_plot(res, y="PIP", b=b)
```
![1762916429195](image/README/1762916429195.png)
