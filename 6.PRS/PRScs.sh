# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *# 
# # # # PRScsx
# 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *# 

# for AJHG PRScs
eur2="/public/home/shilulu/Wulab/sll/ARHL/PRScsx/gwas/ajhg_eur_hm3_prscsx"
ref_dir="/public/share/wchirdzhq2022/Wulab_share/PRScsx/ldblk_1kg_eur"
tag_bim="/public/home/shilulu/Wulab/UKB_data/UKB_Array_genotype/ukb22418_chrall_b0_v2_maf001_hwd1e-10"
gwas_sum="$eur2"
gwas_num="331123"
out_dir="/public/home/shilulu/Wulab/sll/ARHL/PRScsx/chr_result/ARHL_AJHG_prscs"
cmd="python /public/home/shilulu/software/PRScs/PRScs.py \
--ref_dir=$ref_dir \
--bim_prefix=$tag_bim \
--sst_file=$gwas_sum \
--n_gwas=$gwas_num \
--chrom={TASK_ID} \
--out_dir=$out_dir"
qsubshcom "$cmd" 20 100G PRScs 90:00:00 "-array=1-22"

cd chr_result
cat ARHL_AJHG_prscs_pst_eff_a1_b0.5_phiauto_chr*.txt > ../ARHL_AJHG_prscs_pst_eff_a1_b0.5_phiauto_chr_all.txt
# run PRS using plink
bfile="/public/share/wchirdzhq2022/Wulab_share/UKB_data/UKB_ancestry/ukb22418_EUR_b0_v2_maf001_hwd1e-10_grm005"
score="ARHL_AJHG_prscs_pst_eff_a1_b0.5_phiauto_chr_all.txt"
cmd="plink --bfile $bfile --score $score 2 4 6 sum --out EUR_AJHG_prscs_chr_all"
qsubshcom "$cmd" 1 100G plink 90:00:00 ""