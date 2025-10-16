# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *# 
# # # # PRScsx
# 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *# 
SNP          A1   A2   BETA      SE


# run PRScsx
eur1="/public/home/shilulu/Wulab/sll/ARHL/PRScsx/gwas/mvp_eur_hm3_prscsx"
eur2="/public/home/shilulu/Wulab/sll/ARHL/PRScsx/gwas/ajhg_eur_hm3_prscsx"
eas="/public/home/shilulu/Wulab/sll/ARHL/PRScsx/gwas/mvp_bbj_eas_hm3_prscsx"
afr="/public/home/shilulu/Wulab/sll/ARHL/PRScsx/gwas/mvp_afr_hm3_prscsx"
amr="/public/home/shilulu/Wulab/sll/ARHL/PRScsx/gwas/mvp_amr_hm3_prscsx"
ref_dir="/public/share/wchirdzhq2022/Wulab_share/PRScsx"
tag_bim="/public/home/shilulu/Wulab/UKB_data/UKB_Array_genotype/ukb22418_chrall_b0_v2_maf001_hwd1e-10"
gwas_sum="$eur1,$eur2,$eas,$afr,$amr"
gwas_num="406300,331123,184699,107791,52348"
gwas_pop="EUR,EUR,EAS,AFR,AMR"
out_dir="/public/home/shilulu/Wulab/sll/ARHL/PRScsx/chr_result"
out_name="ARHL"
cmd="python /public/home/shilulu/software/PRScsx/PRScsx.py \
--ref_dir=$ref_dir \
--bim_prefix=$tag_bim \
--sst_file=$gwas_sum \
--n_gwas=$gwas_num \
--pop=$gwas_pop \
--chrom={TASK_ID} \
--meta=True \
--out_dir=$out_dir \
--out_name=$out_name"
qsubshcom "$cmd" 20 100G PRScsx 90:00:00 "-array=1-22"

# merge all chr result
cd chr_result
cat ARHL_EUR_pst_eff_a1_b0.5_phiauto_chr*.txt > ../ARHL_EUR_pst_eff_a1_b0.5_phiauto_chr_all.txt
cat ARHL_EAS_pst_eff_a1_b0.5_phiauto_chr*.txt > ../ARHL_EAS_pst_eff_a1_b0.5_phiauto_chr_all.txt
cat ARHL_AFR_pst_eff_a1_b0.5_phiauto_chr*.txt > ../ARHL_AFR_pst_eff_a1_b0.5_phiauto_chr_all.txt
cat ARHL_AMR_pst_eff_a1_b0.5_phiauto_chr*.txt > ../ARHL_AMR_pst_eff_a1_b0.5_phiauto_chr_all.txt
cat ARHL_META_pst_eff_a1_b0.5_phiauto_chr*.txt > ../ARHL_META_pst_eff_a1_b0.5_phiauto_chr_all.txt

# run PRS using plink
for i in EUR AFR EAS SAS; do 
    dir="/public/share/wchirdzhq2022/Wulab_share/UKB_data/UKB_ancestry"
    bfile=$dir/ukb22418_${i}_b0_v2_maf001_hwd1e-10_grm005
    score="ARHL_META_pst_eff_a1_b0.5_phiauto_chr_all.txt"
    cmd="plink --bfile $bfile --score $score 2 4 6 sum --out ${i}_prscsx_META"
    qsubshcom "$cmd" 1 100G plink 90:00:00 ""
done  

