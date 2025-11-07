### 1.`popcorn`  trans-ethnic heritability and genetic correlation
```
popcorn 用于估计在不同群体中常见SNP的因果效应量的相关性，同一性状。该方法纳入了研究中所有SNP的关联信息（区别于其他需要事先进行LD-pruning的方法），
而且仅使用GWAS的概括数据，减少了计算负担，并避免了无法获取基因型数据的问题。
```
#### step 1 计算ld
```
popcorn compute -v 2 --bfile1 Popcorn/test/EUR_ref --bfile2 Popcorn/test/EAS_ref Popcorn/test/EUR_EAS_test_ge.cscore

# --gen_effect 计算遗传效应相关，两步都需加上该参数，无则结果为遗传冲击相关
# 软件提供了 eur-eas的
```
#### step 2 计算correlation
```
popcorn fit -v 1 --cfile EUR_EAS_all_gen_imp.cscore --sfile1 $eur.txt.popcornin --sfile2 $eas.txt.popcornin correlation.txt

# summary data col: rsid/SNP a1/A1 a2/A2 af N beta SE 其中af对应a2的freq
# --gen_effect 计算遗传效应相关，两步都需加上该参数
# --sfile1 和 sfile2顺序需与 cscore 中一致
```

### 2.`LDSC` cross trait or trait cross cohort correlation (same ancestry)
```
conda activate ldsc

#! /bin/bash
sums="/public/home/shilulu/software/ldsc/munge_sumstats.py"
ldsc="/public/home/shilulu/software/ldsc/ldsc.py"
hm3="/public/home/shilulu/Wulab/LDSC/eur_w_ld_chr/w_hm3.snplist"
weight="/public/home/shilulu/Wulab/LDSC/eur_w_ld_chr/"
gwas="gwas.txt"
outdir="/ldsc"

prefix=$(basename $gwas ".txt")

#- munge sumstats
cmd="$sums --sumstats $gwas --merge-alleles $hm3 --chunksize 500000 --a1 A1 --a2 A2 --out $outdir/$prefix"
qsubshcom  "$cmd" 1 10G ldsc_sum 1:00:00 ""

#- correlation
order=(ARHL_meta Tiredness Neuroticism LongStandIllness Overall_health_rating Noisy_workplace Snoring Insomnia Taking_other_prescription_medications Waist_circumference Loneliness derpess_mood Past_tobacco_smoking)
files=""
for trait in "${order[@]}"; do
    res=$trait.sumstats.gz
    files+="$res,"
done
files=${files%,}

cmd="$ldsc --rg $files --ref-ld-chr $weight --w-ld-chr $weight --out Genetic_corr"
qsubshcom "$cmd" 1 20G ldsc_corre 60:00:00 ""
```

### 3.`GSMR` putative causal association between two phenotypes
```
#! /bin/bash 
gcta="/public/home/wuyang/bin/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
ref="/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eur/g1000_eur"
exposure="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/06.trait_association/GSMR/exposure.txt"
outcome="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/06.trait_association/GSMR/outcome.txt"

cmd="$gcta --bfile $ref \
--gsmr-file $exposure $outcome \
--gsmr-direction 2 \
--gwas-thresh 5e-08 \
--diff-freq 0.5 \
--clump-r2 0.05 \
--gsmr2-beta \
--gsmr-snp-min 10 \
--heidi-thresh 0.01 \
--effect-plot \
--out ARHL_gsmr \
--thread-num 20"
qsubshcom "$cmd" 20 100G gsmr 5:00:00 ""
```