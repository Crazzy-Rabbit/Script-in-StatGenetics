run_logistic <- function(pheno, ancestry, profile){
    pc_dir <- "/public/home/shilulu/Wulab/UKB_data/UKB_ancestry/PCA/"
    library(data.table)
    pc = fread(paste0(pc_dir, "ukb22418_", ancestry, "_b0_v2_maf001_hwd1e-10_grm005_pca20.proj.eigenvec"), 
            select=c(1:12), col.names=c("FID", "IID", paste0("PC", 1:10)))
    pc$FID = as.character(pc$FID)
    pc$IID = as.character(pc$IID)
    
    phe = fread(pheno, select = c("FID", "IID", "age", "sex", "has_HL_icd10"))
    phe$FID = as.character(phe$FID)
    phe$IID = as.character(phe$IID)
    phe = phe[phe$age >= 50, ]

    pheno_covar = merge(phe, pc, by=c("FID", "IID"))
    dt_model = pheno_covar[, !c("FID", "IID"), with=FALSE]

    null.model <- glm(has_HL_icd10 ~., data=dt_model, family=binomial)
    null.r2 <- pscl::pR2(null.model)["McFadden"]

    prs.result <- NULL
    p.threshold <- c("5e-8", "1e-5", "0.001", "0.01", "0.05", "0.1", "0.5", "1")

    for(i in p.threshold){
        prs <- fread(paste0(profile, i,".profile"))
        prs$FID = as.character(prs$FID)
        prs$IID = as.character(prs$IID)
        # Merge the prs with the phenotype matrix
        pheno.prs <- merge(pheno_covar, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
        dt_model_prs <- pheno.prs[, !c("FID", "IID"), with=FALSE]

        # perform a logistic regression on trait with PRS and the covariates
        model <- glm(has_HL_icd10 ~., data=dt_model_prs, family=binomial)
        model.r2 <- pscl::pR2(model)["McFadden"]

        prs.r2 <- model.r2-null.r2
        # obtain the coeffcient and p-value of association of PRS as follow
        prs.coef <- summary(model)$coefficients["SCORE",]
        prs.beta <- as.numeric(prs.coef[1])
        prs.se <- as.numeric(prs.coef[2])
        prs.p <- as.numeric(prs.coef[4])

        prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
    }
    return(prs.result)
}


#==========================================================
setwd("/public/home/shilulu/Wulab_project/ARHL/P+T")
pheno <- "/public/home/shilulu/Wulab_project/ARHL/HL_pheno_unrelate_hlAgeSex.pheno"

meta_result <- run_logistic(pheno, "AFR", "AFR/Meta_hm3_")