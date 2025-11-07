#============================================================#
#
#  scRNA expression treatment of Eshel et al
#
#============================================================#
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)

# dir for rds file
scDT <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT"
obj <- readRDS(file.path(scDT, "Eshel_et_al_GSE234926_P90_Seurat.rds"))

#- cell ID to cell type
obj.meta <- obj@meta.data %>% select(celltype)
obj.meta$celltype <- gsub(" ","_", obj.meta$celltype)
obj.meta$cell <- rownames(obj.meta)

#- count
obj.counts <- obj@assays$RNA@counts %>% as.data.frame()
obj.counts$gene <- rownames(obj.counts)
obj.counts.long <- obj.counts %>%
  gather(cell,exp,-gene) %>% 
  left_join(obj.meta, by="cell") %>%
  group_by(gene, celltype) %>%
  summarise(exp=sum(exp)) %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp))

dat <- obj.counts.long
#- keep only mouse to human 1to1 mapped orthologs
dir <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing s/NC_revision"
m2h <- fread(file.path(dir, "mouse_human_homologs.txt"))
names(m2h) <- c("MumSYM", "HumSYM")  

dat <- merge(dat, m2h, by.x="gene", by.y="MumSYM")

#- fill in NAs with 0:
tmp <- dat %>% pivot_wider(id_cols=c(gene, HumSYM), names_from=celltype, values_from=exp, values_fill=0, values_fn=list(exp=sum))
dat <- tmp %>% pivot_longer(cols=-c(gene, HumSYM), names_to="celltype", values_to="exp")

#- remove duplicated genes
tmp.dup <- dat %>% add_count(gene) 

#- remove genes not expressed in any cell type
tmp.noexp <- dat %>% 
  group_by(gene) %>% 
  summarise(sum_exp=sum(exp)) %>%
  filter(sum_exp==0)

dat <- dat %>% filter(!gene %in% tmp.noexp$gene)

#- add up count for duplicated cell types
dat <- dat %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp)) %>%
  ungroup()

# 3. Calculate specificity  
dat <- dat %>%
  group_by(gene) %>%
  mutate(specificity=exp_tpm/sum(exp_tpm)) %>%
  ungroup()

## Keep only genes tested in MAGMA  
# 4. Write MAGMA and LDSC input  
#  Get number of genes that represent 10% of the dataset
n_genes <- length(unique(dat$HumSYM))
n_genes_to_keep <- (n_genes * 0.1) %>% round()

### Get LDSC input top 10%
write_group  = function(df,Cell_type) {
  df <- dplyr::select(df, dplyr::all_of(Cell_type), HumSYM)
  
  dir.create("LDSC", showWarnings=FALSE)
  df %>% dplyr::select(HumSYM) %>% readr::write_tsv(paste0("LDSC/",make.names(unique(df[1])),".bed"), col_names=F)
  invisible(df)
}

ldsc_bedfile <- function(d,Cell_type){
  d_spe <- d %>% group_by(.data[[Cell_type]]) %>% slice_max(order_by=specificity, n=n_genes_to_keep, with_ties=TRUE)
  d_spe %>% group_split(.keep = TRUE) %>% purrr::walk(~ write_group(.x, Cell_type))
}

### Write MAGMA/LDSC input files 
# Filter out genes with expression below 1 TPM.

# change as your out dir
setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Eshel")
dat %>% filter(exp_tpm>1) %>% ldsc_bedfile("celltype")
control <- as.data.frame(unique(dat$HumSYM))
readr::write_tsv(control,  "LDSC/control.bed", col_names=F)





### scDRS 
library(Seurat)
library(stringr)
library(SeuratDisk) 

scDT <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT"
obj <- readRDS(file.path(scDT, "Eshel_et_al_GSE234926_P90_Seurat.rds"))

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method="LogNormalize", verbose=FALSE)

setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT")
SaveH5Seurat(obj, filename = "Eshel_et_al.h5seurat", overwrite = TRUE)
Convert("Eshel_et_al.h5seurat", dest = "h5ad", assay = "RNA", overwrite = TRUE)