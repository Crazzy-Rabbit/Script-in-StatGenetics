# genetrate own annotate file and run ldsc

# # # # # # # # # # # # # # # 
## Step 1: Creating an annot file
# # # 1) generate matrix in python
python -c ""
import pandas as pd 
import scanpy as sc 

SC = sc.read_h5ad("ScRNA-seq-P8_P12_P20_mouse_cochlea_Jun27.h5ad")
expr_df = SC.to_df()
expr_df['cell_type'] = SC.obs['cell_type'].values

# per gene avrg expression in per cell type, rm expr 0 in all cell type
mean_expr_in_each_cell_type = expr_df.groupby('cell_type').mean().T
mean_expr_in_each_cell_type = mean_expr_in_each_cell_type.loc[(mean_expr_in_each_cell_type != 0).any(axis=1)]
mean_expr_in_each_cell_type.to_csv('Sc_P8_12_20.txt', sep="\t")
#================================================================

# # # 2) generate homologous gene in R
# if (!require("BiocManager", quietly = TRUE)) 
#     install.packages("BiocManager")
# BiocManager::install("biomaRt")
library(dplyr)
library(data.table)
library(biomaRt)
listEnsembl()
# connect Ensembl
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
## 构建转化函数
homologs <- getBM(attributes = c("chromosome_name","ensembl_gene_id", "external_gene_name",
                                 "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type"),
                  filters = "with_mmusculus_homolog", 
                  values = TRUE, 
                  mart = human)
# head(homologs)
homologs_1to1 <- subset(homologs, mmusculus_homolog_orthology_type == "ortholog_one2one")
# rm MT X Y chromosome genes
homGene_1to1_chr1_22 <- homologs_1to1 %>% 
  filter(!chromosome_name %in% c("MT", "X", "Y")) %>%
  rename(human_ensembl_gene = ensembl_gene_id)

fwrite(homGene_1to1_chr1_22, file="HomGene1to1_mouse2human.txt", sep="\t", row.names=FALSE)
#================================================================

# # # 3) extract 1:1 homologs scRNA genes from mouse cochlea
# 同时生成control.GeneSet文件，也就是所有的基因
library(dplyr)
library(data.table)
setwd("/public/share/wchirdzhq2022/Wulab_share/LDSC/Mouse_cochleae/")
scGene = fread("Sc_P8_12_20.txt")
homologs_1to1 = fread("HomGene1to1_mouse2human.txt")[,2:3,with=FALSE]
colnames(homologs_1to1) = c("human_gene", "gene")

merge_df = merge(scGene, homologs_1to1, by="gene")
merge_df = subset(merge_df, select = -c(gene))
colnames(merge_df) = gsub("[[:space:]/-]", "_", colnames(merge_df))
colnames(merge_df) = gsub("'", "_", colnames(merge_df))
merge_df = merge_df %>% select(human_gene, everything())

control = merge_df[, 1]
fwrite(control, file="Sc_P8_12_20_control.GeneSet", sep="\t", col.names=FALSE)

lapply(2:ncol(merge_df), function(i){
  temp_dt = data.frame(merge_df[[1]], merge_df[[i]])
  # colnames(df1) = names(merge_df)[c(37, i)]
  sorted_dt = temp_dt[order(temp_dt[[2]], decreasing = TRUE), ]
  top_10 = sorted_dt[1:floor(0.1*nrow(sorted_dt)), ]
  top_10_gene = top_10[, 1, drop=FALSE]

  cell_type = colnames(merge_df)[i]
  fwrite(top_10_gene, file=paste0("Sc_P8_12_20_", cell_type, ".GeneSet"), sep="\t", col.names=FALSE)
})
#================================================================

# cd /public/share/wchirdzhq2022/Wulab_share/LDSC/1000G_EUR_Phase3_plink
# for chr in {1..22}; do
#   echo "1000G.EUR.QC.${chr}" >> merge_list.txt
# done
# plink --bfile 1000G.EUR.QC.1 --merge-list merge_list.txt --make-bed --out ../1000G.EUR.QC
# # # 4) use make_annot.py to generate annot.gz file 
# using 100K widow as Hilary K. Finucane et al. Nat Genet

# generate ENSG.coord.txt
library(biomaRt)
library(data.table)
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
human_gene <- getBM(attributes = c("chromosome_name","ensembl_gene_id", "external_gene_name"),
                    filters = "", values = TRUE, mart = human)
colnames(human_gene) = c("CHR", "ENSG_ID", "GENE")

glist_hg19 <- fread("/public/home/shilulu/script/plot_smr/glist-hg19")
colnames(glist_hg19) = c("CHR", "START", "END", "GENE")
merge_df <- merge(glist_hg19, human_gene, by="GENE")
glist_hg19_ENSG <- merge_df[, c("ENSG_ID", "CHR.y", "START", "END")]
colnames(glist_hg19_ENSG) = c("GENE", "CHR", "START", "END")

fwrite(glist_hg19_ENSG, file="ENSG_coord.txt", sep="\t")

# run make_annot.py 为control 和 每个cell type生成
conda activate ldsc
bfile="/public/share/wchirdzhq2022/Wulab_share/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC"
ldsc="/public/home/shilulu/software/ldsc"

ls *.GeneSet | while read id; do
  annot=$(basename -- "${id}" ".GeneSet")
  # for chr in {1..22}; do 
  cmd="python ${ldsc}/make_annot.py \
      --gene-set-file ${annot}.GeneSet \
      --gene-coord-file ${ldsc}/ENSG_coord.txt \
      --bimfile ${bfile}.{TASK_ID}.bim \
      --windowsize 100000 \
      --annot-file ${annot}.{TASK_ID}.annot.gz"
  # done 
  qsubshcom "$cmd" 1 10G ldsc_anot 1:00:00 "-array=1-22"
done

## Step 2: Computing LD scores with an annot file
bfile="/public/share/wchirdzhq2022/Wulab_share/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC"
ldsc="/public/home/shilulu/software/ldsc"
ldsc_dir="/public/share/wchirdzhq2022/Wulab_share/LDSC"

awk '{if ($1!="SNP") {print $1} }' ${ldsc_dir}/w_hm3.snplist > ${ldsc_dir}/listHM3.txt

ls *.GeneSet | while read id; do
  annot=$(basename -- "${id}" ".GeneSet")
  cmd="python ${ldsc}/ldsc.py --l2 \
    --bfile ${bfile}.{TASK_ID} \
    --print-snps ${ldsc_dir}/listHM3.txt \
    --ld-wind-cm 1 \
    --annot ${annot}.{TASK_ID}.annot.gz \
    --thin-annot \
    --out ${annot}.{TASK_ID}"
  qsubshcom "$cmd" 1 10G ldsc_l2 1:00:00 "-array=1-22"
dones

# generate .ldcts file for my sc data
ls ${dir}/*GeneSet | while read id; do
  tissue=$(basename -- ${id} | sed 's/^Sc_P8_12_20_//;s/.GeneSet$//')
  cell=${dir}/$(basename -- ${id} "GeneSet")
  control=${dir}/"Sc_P8_12_20_control."
  echo "${tissue} ${cell},${control}" >> Sc_P8_12_20_Mouse_cochleae_gene_expr.ldcts
done
