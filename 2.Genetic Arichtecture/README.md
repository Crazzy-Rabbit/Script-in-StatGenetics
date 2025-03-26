### GCTB_SbayesS
```
# We introduced a method (called BayesS) to infer the action of natural selection 
# on the genetic variants underlying a complex trait. 
# stimating the relationship between the variance of SNP effects and MAF 
# Zeng et al. Nat Gent
```
#### 1 step: imputed the summary data
```
#! /bin/bash
gctb="/public/home/shilulu/software/gctb"
ldref="/GCTB_ref/ukbEUR_Imputed"
gwas="gwas.txt"
out_block="/gctb/impu_block/gwas_gctb_impu"

# step1: impute summary data --block 分blcok跑，总共是591个block，效率更高
for i in {1..591}; do
    cmd1="${gctb} --ldm-eigen ${ldref} --gwas-summary ${gwas} --impute-summary --block {i} --out ${i} --thread 10 > ${out_block}{i}.log 2>&1"
    nohup ${i} &
done

out_impu="gwas_gctb.imputed"
# step2: combined the all block file
cd ./impu_block
gctb --gwas-summary gwas_gctb_impu --merge-block-gwas-summary --out ${out_impu} > ../${out_impu}.log 2>&1

mv ${out_impu}.ma ../${out_impu}.ma
```

#### 2 step: run SbayesS
```
# # # 2.8 M snp + baseline annot
ldref="/GCTB_ref/ukb_50k_bigset_2.8M"
outdir="/gctb/bayesS_annot_chr"
annofile="/GCTB_ref/annot_baseline2.2.txt"
gwas_impu="gwas_gctb.imputed.ma"
outprx=$(basename -- ${gwas_impu} ".imputed.ma")

for i in {1..22}; do
  cmd="gctb --sbayes S --ldm ${ldref}/ukb50k_shrunk_chr{TASK_ID}_mafpt01.ldm.sparse --annot ${annofile} --chr {TASK_ID} \
  --gwas-summary ${gwas_impu} --no-mcmc-bin --thread 10 \
  --out ${outdir}/${outprx}_SbayesS_chr{TASK_ID} > ${outdir}/${outprx}_SbayesS_chr{TASK_ID}.log 2>&1"
  nohup ${i} &
done
```