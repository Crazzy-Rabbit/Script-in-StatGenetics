install.packages("pROC")
library(pROC)


run_AUC <- function(pheno, PRS, cov_only=FALSE){
  pred_auto = fread(PRS, select=c(1,2,7), col.names=c("FID", "IID", "PRS"))
  pred_auto[, FID := as.character(FID)]
  pred_auto[, IID := as.character(IID)]
  auto.y = pheno[pred_auto, on=c("FID", "IID"), nomatch=0]
  auto.y[, PRS_std := scale(PRS)] # normalization
  
  if (!cov_only){
    model = paste("PC", 1:10, sep="", collapse="+") %>%
      paste0("has_HL_icd10 ~ PRS_std + sex + age +", .) %>%
      as.formula %>% 
      glm(., data=auto.y, family=binomial())
  } else {
    model = paste("PC", 1:10, sep="", collapse="+") %>%
      paste0("has_HL_icd10 ~ sex + age +", .) %>%
      as.formula %>% 
      glm(., data=auto.y, family=binomial())
  }

  auc <- roc(auto.y$has_HL_icd10, predict(model, type = "response"))
  
  return(auc)
}
setwd("E:/Shi/ALL_of_my_Job/24-28川大华西/2_project_hearing loss/process/PRS/PRS_UKB/demeta/rawmeta")
pheno = fread("14May2025/eur_age50_pheno.txt")
pheno[, FID := as.character(FID)]
pheno[, IID := as.character(IID)]

ajhg="14May2025/LDpred2_ajhg_auto_plot1.txt"
meta="14May2025/LDpred2_meta_auto_plot1.txt"

(meta_auc = run_AUC(pheno, meta, cov_only = FALSE))
ci.auc(meta_auc)

(ajhg_auc = run_AUC(pheno, ajhg, cov_only = FALSE))
ci.auc(ajhg_auc)

(cov_only = run_AUC(pheno, meta, cov_only=TRUE))
ci.auc(cov_only)

roc.test(meta_auc, ajhg_auc, method = "delong")
roc.test(meta_auc, cov_only, method = "delong")
roc.test(ajhg_auc, cov_only, method = "delong")

df_roc <- bind_rows(
  data.frame(
    FPR = rev(cov_only$specificities),
    TPR = rev(cov_only$sensitivities),
    Model = "Covariates-only",
    AUC = as.numeric(cov_only$auc)
  ),
  data.frame(
    FPR = rev(ajhg_auc$specificities),
    TPR = rev(ajhg_auc$sensitivities),
    Model = "PRS2022",
    AUC = as.numeric(ajhg_auc$auc)
  ),
  data.frame(
    FPR = rev(meta_auc$specificities),
    TPR = rev(meta_auc$sensitivities),
    Model = "PRSMeta",
    AUC = as.numeric(meta_auc$auc)
  )
)


p = ggplot(df_roc, aes(x = 1 - FPR, y = TPR, color = Model)) +
  geom_line(linewidth=0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_classic() +
  labs(title = "ROC Curves for PRS Models", x = "False Positive Rate", y = "True Positive Rate") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values = c(
    "Covariates-only" = "#f3dbc7",
    "PRS2022" = "#fc5151",
    "PRSMeta" = "#5c5ca6"
  ))
p

ggsave("AUC.png", p, width=10, height=10, dpi=600)
