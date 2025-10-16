#-----------------------------------------------------------------//
#
# LDpred2 auto model
#
#-----------------------------------------------------------------//
phe <- "/public/home/shilulu/Wulab_project/ARHL/HL_pheno_unrelate_hlAgeSex.pheno"
pc_dir <- "/public/home/shilulu/Wulab/UKB_data/UKB_ancestry/PCA/"

library(bigsnpr)
library(data.table)
library(magrittr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

#------------------ 1. read pheno and covariate file ----//
ukb_phe = fread(phe, select = c("FID", "IID", "age", "sex", "has_HL_icd10"))
ukb_phe[, FID := as.character(FID)]
ukb_phe[, IID := as.character(IID)]
eur = fread("ukb_eur_age50.sample", col.names=c("FID","IID"))
eur[, FID := as.character(FID)]
eur[, IID := as.character(IID)]
eur_phe = merge(eur, ukb_phe, by=c("FID", "IID"))

pc_info = paste0(pc_dir, "ukb22418_EUR_b0_v2_maf001_hwd1e-10_grm005_pca20.proj.eigenvec")
pc = fread(pc_info, select=c(1:12), col.names=c("FID", "IID", paste0("PC", 1:10)))
pc[, FID := as.character(FID)]
pc[, IID := as.character(IID)]

pheno = merge(eur_phe, pc, by=c("FID","IID"))

#------------------ 2. load and transform the summary statisfic file ----//
sums = "/public/home/shilulu/Wulab_project/ARHL/Demeta_raw_hm3_AJHG_rescale_keep_metaonly_SNP.txt"
tmp_gwas = bigreadr::fread2(sums)
names(tmp_gwas) = c("rsid", "a1", "a0", "freq", "beta", "beta_se", "p", "n_eff")

#------------------ 3. calculate LD matrix ----//
prefix = "/public/home/shilulu/Wulab_project/ARHL/LDpred2/ukb_eur_age50_maf.001912"
rds_file = paste0(prefix, ".rds"); 
bed_file = paste0(prefix, ".bed")
if (!file.exists(rds_file)) {
    rds_file = snp_readBed(bed_file)
}
obj.bigSNP = snp_attach(rds_file)

# extract SNP info from genotype
map = obj.bigSNP$map[-3]
names(map) = c("chr", "rsid", "pos", "a1", "a0")
gwas = merge(tmp_gwas, map, by=c("rsid"), suffixes=c(".gwas", ".map"))
gwas = gwas[, c("chr", "pos", "rsid", "a1.gwas", "a0.gwas", "freq", "beta", "beta_se", "p", "n_eff")]
names(gwas) = c("chr", "pos", "rsid", "a1", "a0", "freq", "beta", "beta_se", "p", "n_eff")

info_snp = snp_match(gwas, map)

genotype = obj.bigSNP$genotypes
genotype = snp_fastImputeSimple(genotype)
# CM info from 1000 Genome
CHR = map$chr; 
POS = map$pos; 
SNP_ID = map$rsid 
POS2 = snp_asGeneticPos(CHR, POS, dir="/public/home/shilulu/Wulab/LDpred2")

# LD corr
corr = NULL; 
ld = NULL; 
fam.order = NULL
NCORES = nb_cores()
tmp = tempfile(tmpdir = "/public/home/shilulu/Wulab_project/ARHL/LDpred2/tmp")
on.exit(file.remove(paste0(tmp, ".sbk")), add=TRUE)
for (chr in 1:22){
    ind.chr = which(info_snp$chr == chr)
    ind.chr2 = info_snp$`_NUM_ID_`[ind.chr]

    corr0 = snp_cor(genotype, ind.col=ind.chr2, infos.pos=POS2[ind.chr2], size=3 / 1000, ncores=NCORES)
    if (chr == 1){
        ld = Matrix::colSums(corr0^2)
        corr = as_SFBM(corr0, tmp)
    } else {
        ld = c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}

fam.order = as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))
fam.order[, FID := as.character(FID)]
fam.order[, IID := as.character(IID)]

#------------------ 4. LDSC ----//
df_beta = info_snp[, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
df_beta$chi2 = (df_beta$beta / df_beta$beta_se)^2

ldsc = snp_ldsc(ld, length(ld), chi2=df_beta$chi2, sample_size=df_beta$n_eff, blocks=NULL)
h2_est = ldsc[["h2"]]

#------------------ 5. auto model ----//
NCORES = 4
set.seed(1)
multi_auto = snp_ldpred2_auto(corr,df_beta,h2_init=h2_est,ncores=NCORES,
                                vec_p_init=seq_log(1e-4,0.5, length.out=30),
                                allow_jump_sign=FALSE, shrink_corr=0.95)

#------------------ 6. obtain model PRS ----//
(range = sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep = which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

beta_auto = rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
ind.test = 1:nrow(genotype)
pred_auto = big_prodVec(genotype, beta_auto, ind.row=ind.test, ind.col=info_snp[["_NUM_ID_"]])

auto_plot = data.table(sample=fam.order, prs=pred_auto)
fwrite(auto_plot, "LDpred2_ajhg_auto_plot.txt", sep="\t")

#------------------ 7. calculate R2 ----//
y = pheno[fam.order, on = c("FID", "IID"), nomatch = 0]
null.model = paste("PC", 1:10, sep="", collapse="+") %>%
    paste0("has_HL_icd10 ~ sex + age +", .) %>%
    as.formula %>%
    glm(., data=y, family=binomial())
null.r2 = pscl::pR2(null.model)["McFadden"]


reg.formula = paste("PC", 1:10, sep="", collapse="+") %>%
    paste0("has_HL_icd10 ~ PRS + sex + age +", .) %>%
    as.formula
reg.dat = y
reg.dat$PRS = pred_auto
auto.model = glm(reg.formula, data=reg.dat, family=binomial())
auto.r2 = pscl::pR2(auto.model)["McFadden"]
result = data.table(auto = auto.r2 - null.r2, null = null.r2)

fwrite(result, "LDpred2_ajhg_auto.txt", sep="\t")
#---------------------------------
qsubshcom "Rscript LDpred2_auto_ajhg.r" 20 100G ldpred2 10:00:00 ""

#-----------------------------------------------------------------//
#
# LDpred2 infinitesimal model
#
#-----------------------------------------------------------------//
phe <- "/public/home/shilulu/Wulab_project/ARHL/HL_pheno_unrelate_hlAgeSex.pheno"
pc_dir <- "/public/home/shilulu/Wulab/UKB_data/UKB_ancestry/PCA/"

library(bigsnpr)
library(data.table)
library(magrittr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

#------------------ 1. read pheno and covariate file ----//
ukb_phe = fread(phe, select = c("FID", "IID", "age", "sex", "has_HL_icd10"))
ukb_phe[, FID := as.character(FID)]
ukb_phe[, IID := as.character(IID)]
eur = fread("ukb_eur_age50.sample", col.names=c("FID","IID"))
eur[, FID := as.character(FID)]
eur[, IID := as.character(IID)]
eur_phe = merge(eur, ukb_phe, by=c("FID", "IID"))

pc_info = paste0(pc_dir, "ukb22418_EUR_b0_v2_maf001_hwd1e-10_grm005_pca20.proj.eigenvec")
pc = fread(pc_info, select=c(1:12), col.names=c("FID", "IID", paste0("PC", 1:10)))
pc[, FID := as.character(FID)]
pc[, IID := as.character(IID)]

pheno = merge(eur_phe, pc, by=c("FID","IID"))

#------------------ 2. load and transform the summary statisfic file ----//
sums = "/public/home/shilulu/Wulab_project/ARHL/Demeta_raw_hm3_Meta_rescale_keep_metaonly_SNP.txt"
tmp_gwas = bigreadr::fread2(sums)
names(tmp_gwas) = c("rsid", "a1", "a0", "freq", "beta", "beta_se", "p", "n_eff")

#------------------ 3. calculate LD matrix ----//
prefix = "/public/home/shilulu/Wulab_project/ARHL/LDpred2/ukb_eur_age50_maf.001912"
rds_file = paste0(prefix, ".rds"); 
bed_file = paste0(prefix, ".bed")
if (!file.exists(rds_file)) {
    rds_file = snp_readBed(bed_file)
}
obj.bigSNP = snp_attach(rds_file)

# extract SNP info from genotype
map = obj.bigSNP$map[-3]
names(map) = c("chr", "rsid", "pos", "a1", "a0")
gwas = merge(tmp_gwas, map, by=c("rsid"), suffixes=c(".gwas", ".map"))
gwas = gwas[, c("chr", "pos", "rsid", "a1.gwas", "a0.gwas", "freq", "beta", "beta_se", "p", "n_eff")]
names(gwas) = c("chr", "pos", "rsid", "a1", "a0", "freq", "beta", "beta_se", "p", "n_eff")

info_snp = snp_match(gwas, map)

genotype = obj.bigSNP$genotypes
genotype = snp_fastImputeSimple(genotype)
# CM info from 1000 Genome
CHR = map$chr; 
POS = map$pos; 
SNP_ID = map$rsid 
POS2 = snp_asGeneticPos(CHR, POS, dir="/public/home/shilulu/Wulab/LDpred2")

# LD corr
corr = NULL; 
ld = NULL; 
fam.order = NULL
NCORES = nb_cores()
tmp = tempfile(tmpdir = "/public/home/shilulu/Wulab_project/ARHL/LDpred2/tmp")
on.exit(file.remove(paste0(tmp, ".sbk")), add=TRUE)
for (chr in 1:22){
    ind.chr = which(info_snp$chr == chr)
    ind.chr2 = info_snp$`_NUM_ID_`[ind.chr]

    corr0 = snp_cor(genotype, ind.col=ind.chr2, infos.pos=POS2[ind.chr2], size=3 / 1000, ncores=NCORES)
    if (chr == 1){
        ld = Matrix::colSums(corr0^2)
        corr = as_SFBM(corr0, tmp)
    } else {
        ld = c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}

fam.order = as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))
fam.order[, FID := as.character(FID)]
fam.order[, IID := as.character(IID)]

#------------------ 4. LDSC ----//
df_beta = info_snp[, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
df_beta$chi2 = (df_beta$beta / df_beta$beta_se)^2

ldsc = snp_ldsc(ld, length(ld), chi2=df_beta$chi2, sample_size=df_beta$n_eff, blocks=NULL)
h2_est = ldsc[["h2"]]

#------------------ 5. infinitesimal model ----//
beta_inf = snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
ind.test = 1:nrow(genotype)
pred_inf = big_prodVec(genotype, beta_inf, ind.row=ind.test, ind.col=info_snp[["_NUM_ID_"]])

inf_plot = data.table(sample=fam.order, prs=pred_inf)
fwrite(inf_plot, "LDpred2_meta_inf_plot.txt", sep="\t")

#------------------ 6. calculate R2 ----//
y = pheno[fam.order, on = c("FID", "IID"), nomatch = 0]
null.model = paste("PC", 1:10, sep="", collapse="+") %>%
    paste0("has_HL_icd10 ~ sex + age +", .) %>%
    as.formula %>%
    glm(., data=y, family=binomial())
null.r2 = pscl::pR2(null.model)["McFadden"]


reg.formula = paste("PC", 1:10, sep="", collapse="+") %>%
    paste0("has_HL_icd10 ~ PRS + sex + age +", .) %>%
    as.formula
reg.dat = y
reg.dat$PRS = pred_inf
auto.model = glm(reg.formula, data=reg.dat, family=binomial())
auto.r2 = pscl::pR2(auto.model)["McFadden"]
result = data.table(auto = auto.r2 - null.r2, null=null.r2)

fwrite(result, "LDpred2_meta_inf.txt", sep="\t")
#---------------------------------
qsubshcom "Rscript LDpred2_inf_ajhg.r" 20 100G ldpred2 10:00:00 ""
qsubshcom "Rscript LDpred2_inf_meta.r" 20 100G ldpred2 10:00:00 ""

#-----------------------------------------------------------------//
#
# LDpred2 grid model
#
#-----------------------------------------------------------------//
phe <- "/public/home/shilulu/Wulab_project/ARHL/HL_pheno_unrelate_hlAgeSex.pheno"
pc_dir <- "/public/home/shilulu/Wulab/UKB_data/UKB_ancestry/PCA/"

library(bigsnpr)
library(data.table)
library(magrittr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

#------------------ 1. read pheno and covariate file ----//
ukb_phe = fread(phe, select = c("FID", "IID", "age", "sex", "has_HL_icd10"))
ukb_phe[, FID := as.character(FID)]
ukb_phe[, IID := as.character(IID)]
eur = fread("ukb_eur_age50.sample", col.names=c("FID","IID"))
eur[, FID := as.character(FID)]
eur[, IID := as.character(IID)]
eur_phe = merge(eur, ukb_phe, by=c("FID", "IID"))

pc_info = paste0(pc_dir, "ukb22418_EUR_b0_v2_maf001_hwd1e-10_grm005_pca20.proj.eigenvec")
pc = fread(pc_info, select=c(1:12), col.names=c("FID", "IID", paste0("PC", 1:10)))
pc[, FID := as.character(FID)]
pc[, IID := as.character(IID)]

pheno = merge(eur_phe, pc, by=c("FID","IID"))

#------------------ 2. load and transform the summary statisfic file ----//
sums = "/public/home/shilulu/Wulab_project/ARHL/Demeta_raw_hm3_Meta_rescale_keep_metaonly_SNP.txt"
tmp_gwas = bigreadr::fread2(sums)
names(tmp_gwas) = c("rsid", "a1", "a0", "freq", "beta", "beta_se", "p", "n_eff")

#------------------ 3. calculate LD matrix ----//
prefix = "/public/home/shilulu/Wulab_project/ARHL/LDpred2/ukb_eur_age50_maf.001912"
rds_file = paste0(prefix, ".rds"); 
bed_file = paste0(prefix, ".bed")
if (!file.exists(rds_file)) {
    rds_file = snp_readBed(bed_file)
}
obj.bigSNP = snp_attach(rds_file)

# extract SNP info from genotype
map = obj.bigSNP$map[-3]
names(map) = c("chr", "rsid", "pos", "a1", "a0")
gwas = merge(tmp_gwas, map, by=c("rsid"), suffixes=c(".gwas", ".map"))
gwas = gwas[, c("chr", "pos", "rsid", "a1.gwas", "a0.gwas", "freq", "beta", "beta_se", "p", "n_eff")]
names(gwas) = c("chr", "pos", "rsid", "a1", "a0", "freq", "beta", "beta_se", "p", "n_eff")

info_snp = snp_match(gwas, map)

genotype = obj.bigSNP$genotypes
genotype = snp_fastImputeSimple(genotype)
# CM info from 1000 Genome
CHR = map$chr; 
POS = map$pos; 
SNP_ID = map$rsid 
POS2 = snp_asGeneticPos(CHR, POS, dir="/public/home/shilulu/Wulab/LDpred2")

# LD corr
corr = NULL; 
ld = NULL; 
fam.order = NULL
NCORES = nb_cores()
tmp = tempfile(tmpdir = "/public/home/shilulu/Wulab_project/ARHL/LDpred2/tmp")
on.exit(file.remove(paste0(tmp, ".sbk")), add=TRUE)
for (chr in 1:22){
    ind.chr = which(info_snp$chr == chr)
    ind.chr2 = info_snp$`_NUM_ID_`[ind.chr]

    corr0 = snp_cor(genotype, ind.col=ind.chr2, infos.pos=POS2[ind.chr2], size=3 / 1000, ncores=NCORES)
    if (chr == 1){
        ld = Matrix::colSums(corr0^2)
        corr = as_SFBM(corr0, tmp)
    } else {
        ld = c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}

fam.order = as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))
fam.order[, FID := as.character(FID)]
fam.order[, IID := as.character(IID)]

#------------------ 4. LDSC ----//
df_beta = info_snp[, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
df_beta$chi2 = (df_beta$beta / df_beta$beta_se)^2

ldsc = snp_ldsc(ld, length(ld), chi2=df_beta$chi2, sample_size=df_beta$n_eff, blocks=NULL)
h2_est = ldsc[["h2"]]

#------------------ 5. gred model ----//
(h2_seq = round(ldsc_h2_est * c(0.3, 0.7, 1, 1.4), 4))
(p_seq = signif(seq_log(1e-5, 1, length.out = 21), 2))
(params = expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))

set.seed(1)
NCORES = 4
beta_grid = snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
pred_grid = big_prodMat(genotype, beta_grid, ind.col = info_snp[["_NUM_ID_"]])

library(dplyr)
params %>%
    mutate(sparsity = colMeans(beta_grid == 0), id = row_number()) %>%
    arrange(desc(score)) %>%
    mutate_at(c("score", "sparsity"), round, digits = 3) %>%
    slice(1:10)

best_beta_grid <- params %>%
    mutate(id = row_number()) %>%
    # filter(sparse) %>% 
    arrange(desc(score)) %>%
    slice(1) %>%
    print() %>% 
    pull(id) %>% 
    beta_grid[, .]

ind.test = 1:nrow(genotype)
pred <- big_prodVec(genotype, best_beta_grid, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])


inf_plot = data.table(sample=fam.order, prs=pred)
fwrite(inf_plot, "LDpred2_meta_grid_plot.txt", sep="\t")

#------------------ 6. calculate R2 ----//
y = pheno[fam.order, on = c("FID", "IID"), nomatch = 0]
null.model = paste("PC", 1:10, sep="", collapse="+") %>%
    paste0("has_HL_icd10 ~ sex + age +", .) %>%
    as.formula %>%
    glm(., data=y, family=binomial())
null.r2 = pscl::pR2(null.model)["McFadden"]


reg.formula = paste("PC", 1:10, sep="", collapse="+") %>%
    paste0("has_HL_icd10 ~ PRS + sex + age +", .) %>%
    as.formula
reg.dat = y
reg.dat$PRS = pred_inf
auto.model = glm(reg.formula, data=reg.dat, family=binomial())
auto.r2 = pscl::pR2(auto.model)["McFadden"]
result = data.table(auto = auto.r2 - null.r2, null=null.r2)

fwrite(result, "LDpred2_meta_grid.txt", sep="\t")
#---------------------------------

