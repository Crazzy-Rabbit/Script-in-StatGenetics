# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# employ logistic regression (PRS as x, case/control as y) in each ancestry
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

library(data.table)
library(dplyr)
library(glm2)  # 更稳定的逻辑回归包
library(broom) # extract results of regression

read_file <- function(prs, pheno, pc) {
  # read PRS score
  prs_score <- fread(prs, select = c(1,2,6), col.names = c("FID", "IID", "PRS"))
  prs_score$FID <- as.character(prs_score$FID)
  prs_score$IID <- as.character(prs_score$IID)
  
  # read phe
  ukb_phe <- fread(pheno, select = c("FID", "IID", "age", "sex", "has_HL_icd10"))
  ukb_phe$FID <- as.character(ukb_phe$FID)
  ukb_phe$IID <- as.character(ukb_phe$IID)
  ukb_phe$group <- ukb_phe$has_HL_icd10 # case 1 control 0
  ukb_phe <- ukb_phe[ukb_phe$age >= 50, ]
  
  # read PC info
  pc <- fread(pc, select=c(1:22), col.names=c("FID", "IID", paste0("PC", 1:20)))
  pc$FID <- as.character(pc$FID)
  pc$IID <- as.character(pc$IID)
  
  # combin prs score
  plt_prs <- merge(merge(prs_score, ukb_phe, by=c("FID", "IID")), pc, by=c("FID", "IID"))
  
  return(plt_prs)
}

run_logistic <- function(data){
  # normlization of prs score
  data$PRS_std <- (data$PRS - mean(data$PRS, na.rm = TRUE)) / sd(data$PRS, na.rm = TRUE)
  
  # logistic regression model, covariance are age, sex, top 10 PC
  pc_vars <- paste0("PC", 1:10)
  formula_str <- paste("group ~ PRS_std + age + sex +", paste(pc_vars, collapse = " + "))
  model <- glm(as.formula(formula_str), data=data, family=binomial())
  r2 <- pscl::pR2(model)["McFadden"]
  summary_model <- summary(model)
  
  beta1 <- coef(summary_model)["PRS_std", "Estimate"]
  SE_beta1 <- coef(summary_model)["PRS_std", "Std. Error"]
  
  OR_SD <- exp(beta1)
  CI_lower <- exp(beta1 - 1.96 * SE_beta1)
  CI_upper <- exp(beta1 + 1.96 * SE_beta1)
  
  return(data.frame(OR_SD = OR_SD, CI_lower = CI_lower, CI_upper = CI_upper, R2 = r2))
}



#==========================================================
setwd("/public/share/wchirdzhq2022/Wulab_share/sll/ARHL")

pheno <- "/public/home/shilulu/Wulab_project/ARHL/HL_pheno_unrelate_hlAgeSex.pheno"
pc_dir <- "/public/home/shilulu/Wulab/UKB_data/UKB_ancestry/PCA/"
prs_dir <- "/public/home/shilulu/Wulab_project/ARHL/PRScsx/"


# for (sample in c("EUR", "EAS", "AFR", "SAS")){
for (sample in c("EUR", "EAS", "AFR")){
    pc_info <- paste0(pc_dir, "ukb22418_", sample, "_b0_v2_maf001_hwd1e-10_grm005_pca20.proj.eigenvec")
    # prs <- paste0(prs_dir, sample, "_prscsx_META.profile")
    # prs <- paste0(prs_dir, sample, "_prscsx.profile")
    prs <- paste0(prs_dir, sample, "_AJHG_prscsx.profile")

    data <- read_file(prs, pheno, pc_info)
    result <- run_logistic(data)
    result$ancestry <- sample
    
    # if (!exists("meta_result")){ meta_result <- result } else { meta_result <- rbind(meta_result, result) }
    if (!exists("ajhg_result")){ ajhg_result <- result } else { ajhg_result <- rbind(ajhg_result, result) }
}


meta_result$source <- "meta"
ajhg_result$source <- "ajhg"

# 合并
df <- rbind(meta_result, ajhg_result)
df$ancestry <- factor(df$ancestry, levels = c("EUR", "EAS", "AFR"))

library(ggplot2)

ggplot(df, aes(x = ancestry, y = OR_SD, ymin = CI_lower, ymax = CI_upper, color = source)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(x = "Ancestry", y = "Odds Ratio (per SD increase in PRS)", title = "Comparison of PRS ORs across ancestry and method") +
  scale_color_manual(values = c("meta" = "#1f77b4", "ajhg" = "#ff7f0e"), name = "Analysis") +
  theme_minimal(base_size = 14)+
  geom_text(aes(label = sprintf("%.2f", OR_SD)), 
        position = position_dodge(width = 0.5),
        vjust = -1, size = 3)
# ggsave("forest_plot_meta_vs_ajhg.pdf", width = 8, height = 5)
ggsave("forest_plot_metasingle_vs_ajhg.pdf", width = 8, height = 5)
