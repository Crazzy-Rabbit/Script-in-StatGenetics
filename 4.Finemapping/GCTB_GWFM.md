#### GWFM
###### 1. imputed the summary data
```
#! /bin/bash
gctb="/public/home/shilulu/software/gctb"
ldref="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukbEUR_Imputed"
gwas="/public/home/shilulu/Wulab/sll/ARHL/NC_sup_test/01.meta/All_MVP_Trpchevska_De-Angelis_BBJ_filter.txt"
out_block="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/04.gctb/impu_block/All_MVP_Trpchevska_De-Angelis_BBJ_impu"

# step1: impute summary data --block split block running is a efficiency way, totally 591 block 
qsubshcom "${gctb} --ldm-eigen ${ldref} \
--gwas-summary ${gwas} \
--impute-summary \
--block {TASK_ID} \
--thread 10 \
--out ${out_block} > ${out_block}{TASK_ID}.log 2>&1"  1 100G gctb_impu 1:00:00 "-array=1-591"

out_impu="All_MVP_Trpchevska_De-Angelis_BBJ.imputed"
# step2: combined the all block file to one 
cd ./impu_block
cmd="gctb --gwas-summary All_MVP_Trpchevska_De-Angelis_BBJ_impu --merge-block-gwas-summary --out ${out_impu} > ../${out_impu}.log 2>&1"
qsubshcom "$cmd" 1 100G gctb_impu_merge 1:00:00 ""
```
###### 2. Running GWFM
```
# Run genome-wide fine-mapping analysis
gctb="/public/home/shilulu/software/gctb"
ldref="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukbEUR_Imputed"
annofile="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/annot_baseline2.2.txt"
gwas_impu="All_MVP_Trpchevska_De-Angelis_BBJ.imputed.ma"

outprx=$(basename -- ${gwas_impu} ".imputed.ma")
cmd="${gctb} --sbayes RC --ldm-eigen ${ldref} --annot ${annofile} --gwas-summary ${gwas_impu} --n-dist-auto --write-mcmc-bin --thread 10 \
--out ${outprx}_SbayesRC > ${outprx}_SbayesRC.log 2>&1"
qsubshcom "$cmd" 10 100G sbayesRC 90:00:00 ""

# # # Calculate credible sets
# step1: calculate pairwise ld r2 > 0.5, only need run once
# cmd="gctb --get-ld --ldm-eigen ${ldref} --rsq 0.5 --out LD_rsq05 --thread 10"
# qsubshcom "$cmd" 10 100G gctb_ld 90:00:00 ""
# step2: use ldfile calculate credible sets
gctb="/public/home/shilulu/software/gctb"
ldfile="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukbEUR_Imputed/LD_rsq05.ld.txt"
mcmc_prx="All_MVP_Trpchevska_De-Angelis_BBJ_SbayesRC"

cmd="$gctb --cs \
--ld-file ${ldfile} \
--pip 0.9 \
--mcmc-samples ${mcmc_prx} \
--out ${mcmc_prx}_finemapping > ${mcmc_prx}_finemapping.log 2>&1"
qsubshcom "$cmd" 10 100G sbayesRC 90:00:00 ""
```


#### 可视化结果
###### 单个变异
```
source("/public/home/shilulu/script/plot_GWFM.r")
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL_addchr.gz"
gwfm="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/fine_mapping/ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_SbayesRC.snpRes"
genelist="/public/home/shilulu/script/plot_smr/glist_hg19_refseq.txt"

snp = "rs1126809"

PData = ReadPvalueFromFiles(gwas=gwas, gwfm=gwfm, glist=genelist, windowsize=300000, highlight=snp)
pdf("rs1126809_test.pdf",width = 8,height = 8)
MultiPvalueLocusPlot(data=PData)
dev.off()
```
###### 多个变异循环
```
source("/public/home/shilulu/script/plot_GWFM.r")
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL_addchr.gz"
gwfm="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/fine_mapping/ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_SbayesRC.snpRes"
genelist="/public/home/shilulu/script/plot_smr/glist_hg19_refseq.txt"


SNP = fread("/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/fine_mapping/ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_SbayesRC_finemapping.gcs")
causal = SNP[PIP > 0.9, ][ ,SNP]

causal = data.table(
    rsID = c("rs62623452", "rs13147559", "rs2242416", "rs9493627", "rs35094336", "rs143282422", "rs1126809",
            "rs35887622", "rs2289702", "rs1805005", "rs118174674", "rs61734651", "rs45598239", "rs36062310"),
    CHR = c("2", "4", "6", "6", "8", "10", "11", "13", "15", "16", "18", "20", "21", "22")
)

dir_path = "/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/fine_mapping/plot"
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