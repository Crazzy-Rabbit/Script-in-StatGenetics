#===========================================================================#
#                               *  gsMap  *                                 #
#                             *             *                               #
#                           *                 *                             #
#      Genetically informed spatial mapping of cells for complex traits     #
#                                                                           #
#    Liyang Song, Wenhao Chen, Junren Hou, Minmin Guo, Jian Yang (2024)     #
#  Spatially resolved mapping of cells associated with human complex traits #
#===========================================================================#
conda activate gsMap

# 01. data format
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL"
outdir="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gsMap"
outprx=$(basename -- "$gwas")
gsmap format_sumstats --sumstats ${gwas} --snp SNP --a1 A1 --a2 A2 --frq freq --beta beta --se SE --p p --n N \
--out ${outdir}/${outprx} > ${outdir}/${outprx}.log 2>&1

# # # Quick Mode
# E16.5_E1S1.MOSTA 
# Embryonic-day, Mouse id, ST section id e.g., E10.5_E1S1: Embryonic-day-10, Mouse-1, Section-1
# ST_sample="/public/share/wchirdzhq2022/Wulab_share/gsMap/gsMap_example_data/ST/E16.5_E1S1.MOSTA.h5ad"
# ST_name="E16.5_E1S1.MOSTA"
# gsMap_resource="/public/share/wchirdzhq2022/Wulab_share/gsMap"
# homo="/public/share/wchirdzhq2022/Wulab_share/gsMap/homologs/mouse_human_homologs.txt"
# workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
# gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gsMap/ARHL_MVP_AJHG_BBJ_reformatMETAL.sumstats.gz"
# gsmap quick_mode \
#     --workdir ${workdir} \
#     --homolog_file ${homo} \
#     --sample_name ${ST_name} \
#     --gsMap_resource_dir ${gsMap_resource} \
#     --hdf5_path ${ST_sample} \
#     --annotation 'annotation' \
#     --data_layer 'count' \
#     --sumstats_file ${gwas} \
#     --trait_name ARHL
# # #
# 有报错，这个mode有问题，numpy版本的问题，作者已更新依赖
# ValueError: numpy.dtype size changed, may indicate binary incompatibility. Expected 96 from C header, got 88 from PyObject
# # # 

# # # step by step recommand this 
# 前三步的文件可以保留，以后进行分析的时候不必重新跑，直接用就可以了，从第4步直接开始，
# 但是，请注意，你设定的workdir要和前三步生成文件的那个大目录相同，也就是结果文件也在这个地方
# 感觉没有out参数设定输出就很难受，不过也可能是每一步的输入比较繁琐，统一工作路径就好做
# 并且前三步产生的结果文件有48G，所以，我放公共目录下了
# 因此，切记workdir别改，把后面的trait name改了就行，最后结果会在report文件夹中，找你trait name对应的文件夹
# 1. find latent representations
# The --workdir parameter specifies the working directory for gsMap, where all outputs will be saved.
ST_name="E16.5_E1S1.MOSTA"
ST_sample="/public/share/wchirdzhq2022/Wulab_share/gsMap/gsMap_example_data/ST/E16.5_E1S1.MOSTA.h5ad"
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
gsmap run_find_latent_representations \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --input_hdf5_path ${ST_sample} \
    --annotation 'annotation' \
    --data_layer 'count'

# 2. generate gene specificity scores
# Objective: Identify homogeneous spots for each spot based on their latent representations, 
# and then generate gene specificity scores (GSS) for each spot by aggregating information from its homogeneous spots.
# 识别空间同质性区域，并形成微区
homo="/public/share/wchirdzhq2022/Wulab_share/gsMap/homologs/mouse_human_homologs.txt"
ST_name="E16.5_E1S1.MOSTA"
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
gsmap run_latent_to_gene \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --annotation 'annotation' \
    --latent_representation 'latent_GVAE' \
    --num_neighbour 51 \
    --num_neighbour_spatial 201 \
    --homolog_file ${homo}

# 3. generate ldscore
# Assign gene specificity scores (GSS) to SNPs and compute the stratified LD score.
# Three SNP to gene linking strategies are available:
# 1). Use TSS
# This strategy uses TSS to assign GSS to SNPs.
hm3_SNP="/public/share/wchirdzhq2022/Wulab_share/gsMap/LDSC_resource/hapmap3_snps/hm"
hg37_annot="/public/share/wchirdzhq2022/Wulab_share/gsMap/genome_annotation/gtf/gencode.v39lift37.annotation.gtf"
bfile="/public/share/wchirdzhq2022/Wulab_share/gsMap/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC"
ST_name="E16.5_E1S1.MOSTA"
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
cmd="gsmap run_generate_ldscore \
        --workdir ${workdir} \
        --sample_name ${ST_name} \
        --chrom {TASK_ID} \
        --bfile_root ${bfile} \
        --keep_snp_root ${hm3_SNP} \
        --gtf_annotation_file ${hg37_annot} \
        --gene_window_size 50000"
qsubshcom "$cmd" 10 100G gsmap 90:00:00 "-array=1-22"


### for conditional analysis
# 1. extract spots of lung and generate .feather file

import os 
import scanpy as sc 
import pandas as pd
import numpy as np
from pathlib import Path

os.chdir("/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo/E16.5_E1S1.MOSTA/find_latent_representations/")
adata = sc.read_h5ad("E16.5_E1S1.MOSTA_add_latent.h5ad", backed='r')

os.chdir("/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo/E16.5_E1S1.MOSTA/latent_to_gene")
df = pd.read_feather("E16.5_E1S1.MOSTA_gene_marker_score.feather")
df.index = df.HUMAN_GENE_SYM

df_subset = pd.DataFrame(df[adata.obs_names[adata.obs.annotation == 'Lung']].mean(axis=1), columns=['Lung'])
df_subset.reset_index(inplace=True)
"
df_subset.rename(columns={"index": "HUMAN_GENE_SYM"}, inplace=True)
"
pth = Path("/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo/Lung/latent_to_gene/Lung_gene_marker_score.feather")
pth.parent.mkdir(parents=True, exist_ok=True)
df_subset.to_feather(pth)


# 2. generate ldscore for GSS of lung
hm3_SNP="/public/share/wchirdzhq2022/Wulab_share/gsMap/LDSC_resource/hapmap3_snps/hm"
hg37_annot="/public/share/wchirdzhq2022/Wulab_share/gsMap/genome_annotation/gtf/gencode.v39lift37.annotation.gtf"
bfile="/public/share/wchirdzhq2022/Wulab_share/gsMap/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC"
ST_name="Lung"
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
cmd="gsmap run_generate_ldscore \
        --workdir ${workdir} \
        --sample_name ${ST_name} \
        --chrom {TASK_ID} \
        --bfile_root ${bfile} \
        --keep_snp_root ${hm3_SNP} \
        --gtf_annotation_file ${hg37_annot} \
        --gene_window_size 50000"
qsubshcom "$cmd" 10 100G gsmap 90:00:00 "-array=1-22"

# 3. rename the LD score of lung as additional baseline
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
cd ${workdir}/Lung/generate_ldscore/Lung_chunk1
ls | while read file;
do  
    name=$(echo ${file} | awk -F "Lung." '{print $2}')
    mv ${file} baseline.${name}
done  

# 4: Generate ldscore for each cell and copy the LD score of lung to additional baseline
hm3_SNP="/public/share/wchirdzhq2022/Wulab_share/gsMap/LDSC_resource/hapmap3_snps/hm"
hg37_annot="/public/share/wchirdzhq2022/Wulab_share/gsMap/genome_annotation/gtf/gencode.v39lift37.annotation.gtf"
bfile="/public/share/wchirdzhq2022/Wulab_share/gsMap/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC"
addition_baseline="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo/Lung/additional_baseline"
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
ST_name="E16.5_E1S1.MOSTA"

cmd="gsmap run_generate_ldscore \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --chrom {TASK_ID} \
    --bfile_root ${bfile} \
    --keep_snp_root ${hm3_SNP} \
    --gtf_annotation_file ${hg37_annot} \
    --gene_window_size 50000"
qsubshcom "$cmd" 10 100G run_generate_ldscore 90:00:00 "-array=1-22"

mkdir ${workdir}/E16.5_E1S1.MOSTA/generate_ldscore/additional_baseline
cp ${workdir}/Lung/generate_ldscore/Lung_chunk1/* ${workdir}/E16.5_E1S1.MOSTA/generate_ldscore/additional_baseline/

# # delete additional_baseline folder when complete conditional analysis
# # and you can replicate above cmd in next conditional analysis
rm -r ${workdir}/E16.5_E1S1.MOSTA/generate_ldscore/additional_baseline

# 5: spatial LDSC
# workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
# ST_name="E16.5_E1S1.MOSTA"
workdir="/public/home/shilulu/Wulab/gsMap/Mouse_Embryo/output"
ST_name="E16.5_E1S5.MOSTA"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gsMap/ARHL_MVP_AJHG_BBJ_reformatMETAL.sumstats.gz"
hm3_w="/public/share/wchirdzhq2022/Wulab_share/gsMap/LDSC_resource/weights_hm3_no_hla/weights."
cmd="gsmap run_spatial_ldsc \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL_condition' \
    --sumstats_file ${gwas} \
    --w_file ${hm3_w} \
    --num_processes 30 > spatial_ldsc_condition.log 2>&1"
qsubshcom "$cmd" 10 100G gsmap 10:00:00 ""


# cauchy combination
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
ST_name="E16.5_E1S1.MOSTA"
cmd="gsmap run_cauchy_combination \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL_condition' \
    --annotation 'annotation'"
qsubshcom "$cmd" 1 100G gsmap 10:00:00 ""

# report result
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
ST_name="E16.5_E1S1.MOSTA"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gsMap/ARHL_MVP_AJHG_BBJ_reformatMETAL.sumstats.gz"
cmd="gsmap run_report \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL_condition' \
    --annotation 'annotation' \
    --sumstats_file ${gwas} \
    --top_corr_genes 50 > gsmap_report_result.log 2>&1"
qsubshcom "$cmd" 1 100G gsmap 10:00:00 ""