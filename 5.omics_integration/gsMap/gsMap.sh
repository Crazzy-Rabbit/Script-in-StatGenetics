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
gsmap format_sumstats --sumstats ${gwas} --snp SNP --a1 A1 --a2 A2 --freq freq --beta beta --se SE --p p --n N \
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
# 解释一下，就是如果使用的数据不是人的，那么就得跑这一步找同源基因了，找完在进行下面的操作
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

# 2). Use Enhancer-Gene Linking
cmd="gsmap run_generate_ldscore \
        --workdir ${workdir} \
        --sample_name ${ST_name} \
        --chrom {TASK_ID} \
        --bfile_root ${bfile} \
        --keep_snp_root ${hm3_SNP} \
        --gtf_annotation_file ${hg37_annot} \
        --enhancer_annotation_file 'gsMap_resource/genome_annotation/enhancer/by_tissue/ALL/ABC_roadmap_merged.bed' \
        --snp_multiple_enhancer_strategy 'max_mkscore' \
        --gene_window_enhancer_priority 'enhancer_only'"
qsubshcom "$cmd" 10 100G gsmap 90:00:00 "-array=1-22"

# 3). Use TSS and Enhancer-Gene Linking
# This strategy uses both TSS and enhancer-gene linking to assign GSS to SNPs. 
cmd="gsmap run_generate_ldscore \
        --workdir './example/Mouse_Embryo' \
        --sample_name 'E16.5_E1S1.MOSTA' \
        --chrom {TASK_ID} \
        --bfile_root 'gsMap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC' \
        --keep_snp_root 'gsMap_resource/LDSC_resource/hapmap3_snps/hm' \
        --gtf_annotation_file 'gsMap_resource/genome_annotation/gtf/gencode.v39lift37.annotation.gtf' \
        --gene_window_size 50000 \
        --enhancer_annotation_file 'gsMap_resource/genome_annotation/enhancer/by_tissue/ALL/ABC_roadmap_merged.bed' \
        --snp_multiple_enhancer_strategy 'max_mkscore' \
        --gene_window_enhancer_priority 'gene_window_first'"
qsubshcom "$cmd" 10 100G gsmap 90:00:00 "-array=1-22"

# 4. spatial ldsc
# Objective: Run spatial LDSC to associate spots with traits.
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
ST_name="E16.5_E1S1.MOSTA"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gsMap/ARHL_MVP_AJHG_BBJ_reformatMETAL.sumstats.gz"
hm3_w="/public/share/wchirdzhq2022/Wulab_share/gsMap/LDSC_resource/weights_hm3_no_hla/weights."
gsmap run_spatial_ldsc \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL' \
    --sumstats_file ${gwas} \
    --w_file ${hm3_w} \
    --num_processes 10

# 5. cauchy combination (optional)
# Objective: Aggregate P values of individual spots within specific spatial regions (cell types) to 
# evaluate the association of these regions (cell types) with the trait.
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
ST_name="E16.5_E1S1.MOSTA"
gsmap run_cauchy_combination \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL' \
    --annotation 'annotation'

# 6. report generation (optional)
# Objective: Generate gsMap reports, including visualizations of mapping results and diagnostic plots.
# The default genes for visualization are the top 50 genes 
# whose GSS(gene speicefic score) shows the highest correlation with the -log10 p-values of the trait-cell associations. 
# To select specificity  genes for visualization, use the --selected_genes parameter.
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
ST_name="E16.5_E1S1.MOSTA"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gsMap/ARHL_MVP_AJHG_BBJ_reformatMETAL.sumstats.gz"
gsmap run_report \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL' \
    --annotation 'annotation' \
    --sumstats_file ${gwas} \
    --top_corr_genes 50
