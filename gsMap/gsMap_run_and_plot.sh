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
cmd="gsmap run_report \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL' \
    --annotation 'annotation' \
    --sumstats_file ${gwas} \
    --top_corr_genes 50 \
    --fig_style dark"
qsubshcom "$cmd" 1 80G gsmap 10:00:00 ""


# # # # # plot cauchy combination p value
library(ggplot2)
library(data.table)

data <- fread("gsMap_cauchy_combination_result.txt")
data <- data[order(data$P_Cauchy, decreasing=TRUE),]

# color for dif tissue
custom_colors <- c("Hair cells" = "#DD7694", 
                   "Supporting cells" = "#BCD4E7", 
                   "Surrounding structures" = "#056E83", 
                   "Lateral wall" = "#E9D9BF", 
                   "Circulating cells" = "#D4920A", 
                   "Glial cells" = "#5AA4AE", 
                   "Neurons" = "#65472F")   

# set Name as factor
data$Annotation <- factor(data$Annotation, levels=unique(data$Annotation))

p=ggplot(data, aes(x=Annotation, y=-log10(P_Cauchy), fill=Annotation)) +
  geom_bar(stat="identity", position="dodge", width=0.8, fill="#056E83")+
  scale_y_continuous(expand=c(0,0), limits=c(0, 9))+
  labs(x=NULL, y="-log10(P Cauchy)", title="")+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        legend.position="none",
        legend.title=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        strip.text = element_blank(),
        strip.background=element_blank())
ggsave("plot cauchy combination p value.png", p, width=3.5, height=5, dpi=500)


# # # # # # # # # # # # # # # # # # 
# plot gsMap result for user DIY 
# # # # # # # # # # # # # # # # # # 
library(ggplot2)
library(RColorBrewer)
library(data.table)

data <- fread("E16.5_E1S1.MOSTA_ARHL_gsMap_plot.csv")
my_color <- rev(brewer.pal(5, "RdBu"))

# plot spatial RNA alta
p = ggplot(data, aes(x = sx, y = sy)) +
  geom_point(aes(color = logp)) + 
  scale_color_gradientn(colors = c("#39489f","#39bbec","#F7F7F7", "#F4A582","#b81f25"),
                        limits = c(0, max(data$logp)),
                        breaks = seq(0, max(data$logp), by = 5)) +
  theme_minimal() +
  labs(title = "", x = "", y = "", 
       color = expression(-log[10](p-value))) +
  theme(panel.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank(),
        panel.border = element_blank()) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.position = c(0,0),
        legend.justification = c(0, 0),
        legend.text = element_text(color = "white"),
        legend.title = element_blank(),
        legend.direction = "horizontal") # +
 # facet_wrap(~annotation) 按照组织分类画
ggsave("E16.5_E1S1.MOSTA_ARHL_gsMap_plot.png", p, width=4.5, height=6, dpi=500)

# plot tissue region
my_color <-  c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF1493", 
               "#FF4500", "#228B22", "#8A2BE2", "#FFD700", "#00FFFF", 
               "#FF6347", "#D2691E", "#C71585", "#ADFF2F", "#FF8C00", 
               "#8B0000", "#00BFFF", "#FF7F50", "#4B0082", "#20B2AA", 
               "#B22222", "#F0E68C", "#9932CC", "#2E8B57", "#A52A2A", 
               "#7FFF00")

p = ggplot(data, aes(x = sx, y = sy)) +
  geom_point(aes(color = annotation)) + 
  scale_color_manual(values=my_color) +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(panel.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank(),
        panel.border = element_blank()) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.position = "right",
        legend.text = element_text(color = "black"),
        legend.title = element_blank()) 
ggsave("E16.5_E1S1.MOSTA_ARHL_cell_region.png", p, width=7.5, height=5.6, dpi=500)
