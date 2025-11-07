conda activate ldsc
#*************************#
##       LDSC-SEG        ##
#*************************#
#! /bin/bash
sums="~/software/ldsc/munge_sumstats.py"
ldsc="/public/home/shilulu/software/ldsc"
hm3="~/Wulab/LDSC/eur_w_ld_chr/w_hm3.snplist"
baseline="~/Wulab/LDSC/1000G_EUR_Phase3_baseline/baseline."
weight="~/Wulab/LDSC/weights_hm3_no_hla/weights."
cts="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/08.scEnrich/Eshel/LDSC/config/scRNA_gene_expr.ldcts"
gwas="All_MVP_Trpchevska_De-Angelis_BBJ_filter.gz"

prefix=$(basename $gwas ".gz")

#- munge sumstats
cmd="$sums --sumstats $gwas --merge-alleles $hm3 --chunksize 500000 --a1 A1 --a2 A2 --out $prefix"
qsubshcom "$cmd" 1 100G sumstats 2:00:00 ""

#- run ldsc-seg
cmd="$ldsc --h2-cts ${prefix}.sumstats.gz --ref-ld-chr $baseline --out $prefix --ref-ld-chr-cts $cts --w-ld-chr $weight"
qsubshcom "$cmd" 1 100G ldsc 2:00:00 ""
