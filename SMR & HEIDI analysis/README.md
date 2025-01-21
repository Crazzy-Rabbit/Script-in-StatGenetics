### 1、`smr` & `heidi`分析
##### 1 run SMR and HEIDI
```
smr -bfile bfile --gwas-summary GWAS --beqtl-summary eQTL --out --diff-freq-prop 0.5 --maf 0.01 --cis-wind 2000 --out OUTPUT --thread-num 10
```
`--bfile` SNP基因型数据
`--gwas-summary` summary-level data from GWAS. 格式如下
```
SNP A1 A2 freq b se p N 
rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830

每一列含义如下：
SNP, the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size.
```
##### 2 run generate plot file
```
smr -bfile bfile --gwas-summary GWAS --beqtl-summary eQTL --out --probe cg02853497 --probe-wind 2000 Kb --gene-list Genelist --out OUTPUT --thread-num 10
```
##### plot
```
source("/public/home/shilulu/script/plot_smr/plot_epiSMR.r")
# Read the data file in R:
SMRData = ReadSMRData("chr1.cg02853497.txt")
# Plot the SMR results in a genomic region centred around a probe:
pdf('chr1.cg02853497.pdf',width = 10,height = 8)
SMRLocusPlot(data=SMRData, smr_thresh=5.37e-7, heidi_thresh=0.01, plotWindow=1000, max_anno_probe=8)
dev.off()
```
![image](https://github.com/user-attachments/assets/39be54ef-4f45-4347-b4d9-0af5fbe19640)

### 2、`smr`将多个某QTL的`besd`文件合并
```
ls bl_mqtl_chr*.besd | while read id; do 
    gwas=$(basename -- ${id} ".besd")
    echo `pwd`/"$gwas" >> bl_mqtl_merge_besd.txt
done

smr \
--besd-flist bl_mqtl_merge_besd.txt \
--mecs \
--make-besd-dense \
--thread-num 8 \
--out LBC_BSGS_meta_all
```

### 3、`smr` analysis of two molecular traits
```
smr --bfile mydata --beqtl-summary myexposure --beqtl-summary myoutcome  --out myomics

# --extract-exposure-probe extracts a subset of exposure probes for analysis.
# --extract-outcome-probe extracts a subset of outcome probes for analysis.
# --exclude-exposure-probe excludes a subset of exposure probes from analysis.
# --exclude-outcome-probe excludes a subset of outcome probes from analysis.
```
