# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *# 
#    SBayesRC
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *# 
# /public/home/shilulu/Wulab_project/ARHL/SBayesRC
## 1. impute data using hm3 SNPs
# step1: impute summary data
ldref="/public/home/shilulu/Wulab/GCTB_ref/ukbEUR_HM3"
gwas="Demeta_raw_hm3_AJHG_rescale_keep_metaonly_SNP.txt"
out_block="./impu_block/Demeta_raw_hm3_AJHG_rescale_keep_metaonly_SNP.gctb_impu"

cmd1="gctb --ldm-eigen ${ldref} --gwas-summary ${gwas} --impute-summary --block {TASK_ID} --out ${out_block} --thread 10 > ${out_block}{TASK_ID}.log 2>&1"
qsubshcom "$cmd1" 1 100G gctb_impu 90:00:00 "-array=1-591"

# step2: combined the all block file
out_impu="Demeta_raw_hm3_AJHG_rescale_keep_metaonly_SNP.gctb.imputed"
cd ./impu_block
gctb --gwas-summary Demeta_raw_hm3_AJHG_rescale_keep_metaonly_SNP.gctb_impu --merge-block-gwas-summary --out ${out_impu}
cp ${out_impu}.ma ../${out_impu}.ma

## 2. rerun SbayesRC
ldref="/public/home/shilulu/Wulab/GCTB_ref/ukbEUR_HM3"
annofile="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/annot_baseline2.2.txt"
gwas="Demeta_raw_hm3_AJHG_rescale_keep_metaonly_SNP.gctb.imputed.ma"
outprx=$(basename -- ${gwas} ".imputed.ma")
cmd="gctb --sbayes RC --ldm-eigen ${ldref}  --gwas-summary ${gwas} --annot ${annofile} --thread 20 --out ${outprx}_SbayesRC > ${outprx}_SbayesRC.log 2>&1"
# --write-mcmc-bin  when delete annot use --sbayes R
qsubshcom "$cmd" 10 100G sbayesRC 90:00:00 ""


bfile="/public/share/wchirdzhq2022/Wulab_share/UKB_data/UKB_ancestry/ukb22418_EUR_b0_v2_maf001_hwd1e-10_grm005"
score="Demeta_raw_hm3_AJHG_rescale_keep_metaonly_SNP.gctb_SbayesRC.snpRes"

cmd="plink --bfile $bfile --score $score 2 5 8 header sum --out EUR_AJHG_hm3_SNP"
qsubshcom "$cmd" 1 100G plink 90:00:00 ""