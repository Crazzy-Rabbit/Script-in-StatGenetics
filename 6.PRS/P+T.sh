#!/bin/bash
# run prs
bfile="/public/home/shilulu/Wulab/UKB_data/Array/ukb22418_chrall_b0_v2_maf001_hwd1e-10"
sample="/public/home/shilulu/tmp/PRS_NPDs/pheno"
score="--score PRS_SNP_info_1010.txt 1 4 6 header"

declare -A phenos=(
    [AD]='AD_combined.phe'
    [AN]='AN_combined.phe'
    [BD]='BD_combined.phe'
    [PD]='PD_combined.phe'
    [SCZ]='SCZ_combined.phe'
)

for pheno in "${!phenos[@]}"; do
    file="${phenos[$pheno]}"
    echo "Processing $pheno → $file"
    awk -v OFS="\t" 'NR > 1 {print $1, $2}' ${sample}/${file} > ${sample}/${pheno}.sample

    qsubshcom "plink --bfile $bfile --keep ${sample}/${pheno}.sample $score --threads 30 --out ${pheno}" \
    1 100G plink 5:00:00 ""
done 



# PRS fit
library(data.table)
library(pscl)
library(pROC)

run_logistic <- function(pheno, PRS, covars=c("age", "sex", "assessment_centre")){
    # -------------------------
    # 1️⃣ 读取表型数据
    # -------------------------
    phe = fread(pheno)
    phe <- fread(pheno, select = c("FID", "IID", covars, "has_trait"))
    phe$FID = as.character(phe$FID); phe$IID = as.character(phe$IID)
    phe$has_trait <- as.numeric(as.character(phe$has_trait))

    n_before <- nrow(phe)
    phe <- phe[complete.cases(phe[, c(covars, "has_trait"), with=FALSE]), ]
    cat("Cleaned", n_before - nrow(phe), "rows with missing covariates or has_trait\n")

    # -------------------------
    # 2️⃣ 构建 null model
    # -------------------------
    formula_null <- as.formula(
        paste("has_trait ~", paste(covars, collapse = " + "))
    )
    null.model <- glm(formula_null, data=phe, family=binomial)
    r2.list <- pscl::pR2(null.model)
    null.r2 <- if ("Nagelkerke" %in% names(r2.list)) {
        r2.list[["Nagelkerke"]]
        } else if ("McFadden" %in% names(r2.list)) {
            r2.list[["McFadden"]]
            } else {
                NA
                }

    # -------------------------
    # 3️⃣ 加载 PRS 文件并合并
    # -------------------------
    prs <- fread(PRS, select = c("FID", "IID", "SCORE"))
    prs$FID = as.character(prs$FID); prs$IID = as.character(prs$IID)
    phe_prs = merge(phe, prs[, c("FID", "IID", "SCORE")], by=c("FID", "IID"))
    phe_prs$SCORE = scale(phe_prs$SCORE)

    # -------------------------
    # 4️⃣ 构建包含 PRS 的 model
    # -------------------------
    formula_full <- as.formula(
        paste("has_trait ~ SCORE +", paste(covars, collapse = " + "))
    )
    model = glm(formula_full, data=phe_prs, family=binomial)
    r2.list.full <- pscl::pR2(model)
    model.r2 <- if ("Nagelkerke" %in% names(r2.list.full)) {
        r2.list.full[["Nagelkerke"]]
        } else if ("McFadden" %in% names(r2.list.full)) {
            r2.list.full[["McFadden"]]
            } else {
                NA
                }
    prs.r2 <- model.r2-null.r2

    # -------------------------
    # 5️⃣ 提取 PRS 的效应值
    # -------------------------
    prs.coef <- summary(model)$coefficients["SCORE",]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])

    # ---------------------------
    # 6️⃣ 计算 AUC（预测能力）
    # ---------------------------
    pred.prob <- predict(model, type="response")
    roc.obj <- pROC::roc(phe_prs$has_trait, pred.prob)
    auc.value <- as.numeric(pROC::auc(roc.obj))

    prs.result <- data.frame(PHENO=basename(pheno), P=prs.p, BETA=prs.beta, SE=prs.se, R2=prs.r2, AUC=auc.value)
    
    return(prs.result)
}

# ===================================
# 清洗表型文件
# ===================================
clean_phe <- function(pheno) {
    phe = fread(pheno)
    cat("读取文件:", pheno, "\n")
    cat("原始样本数:", nrow(phe), "\n")

    phe_clean = phe[complete.cases(phe), ]
    
    n_removed = nrow(phe) - nrow(phe_clean)
    cat("已删除某一pheno含 NA 的样本数:", n_removed, "\n")

    # ---------- 统计 case / control ----------
    if (!"has_trait" %in% names(phe_clean)){
        stop("列 'has_trait' 未找到，请检查表型文件！")
    }
    phe_clean$has_trait <- as.numeric(as.character(phe_clean$has_trait))
    n_total <- nrow(phe_clean)
    n_case <- sum(phe_clean$has_trait == 1, na.rm = TRUE)
    n_control <- sum(phe_clean$has_trait == 0, na.rm = TRUE)

    cat("清洗后样本数:", n_total, "\n")
    cat("病例(case=1):", n_case, " (", round(100 * n_case / n_total, 2), "%)\n", sep = "")
    cat("对照(control=0):", n_control, " (", round(100 * n_control / n_total, 2), "%)\n\n", sep = "")

    return(phe_clean)
}

setwd("/public/home/shilulu/tmp/PRS_NPDs")
phe_files <- c(
    "AD_combined.phe",
    "AN_combined.phe",
    "BD_combined.phe",
    "PD_combined.phe",
    "SCZ_combined.phe"
)

for (phe in phe_files) {
    pheno = file.path("pheno", phe)
    phe_qc = clean_phe(pheno)
    fwrite(phe_qc, file=file.path("pheno", paste0(phe, ".delNA")), sep="\t")
}

# ===================================
# 手动输入文件路径（按一一对应顺序）
# ===================================
setwd("/public/home/shilulu/tmp/PRS_NPDs")
phe_files <- c(
    "AD_combined.phe.delNA",
    "AN_combined.phe.delNA",
    "BD_combined.phe.delNA",
    "PD_combined.phe.delNA",
    "SCZ_combined.phe.delNA"
)

prs_files <- c(
    "AD.profile",
    "AN.profile",
    "BD.profile",
    "PD.profile",
    "SCZ.profile"
)

# ===================================
# 批量运行
# ===================================
#- 性别、年龄、研究中心
covars = c("age", "sex", "assessment_centre")
all_results <- data.frame()
for (i in seq_along(phe_files)) {
    phe_file = file.path("pheno", phe_files[i])
    cat("===== Running:", phe_file, "with", prs_files[i], "=====\n")
    res <- run_logistic(pheno = phe_file, PRS = prs_files[i],
                        covars = covars)
    all_results <- rbind(all_results, res)
}

fwrite(all_results, "NPDs_PRS_logistic_age-sex-centre.csv", sep="\t")

#- 性别、年龄、研究中心、BMI
covars = c("age", "sex", "assessment_centre", "BMI")
all_results <- data.frame()
for (i in seq_along(phe_files)) {
    phe_file = file.path("pheno", phe_files[i])
    cat("===== Running:", phe_file, "with", prs_files[i], "=====\n")
    res <- run_logistic(pheno = phe_file, PRS = prs_files[i],
                        covars = covars)
    all_results <- rbind(all_results, res)
}

fwrite(all_results, "NPDs_PRS_logistic_age-sex-centre-BMI.csv", sep="\t")

#- 性别、年龄、研究中心、BMI、饮酒、抽烟、身体活动
covars = c("age", "sex", "assessment_centre", "BMI", "alcohol_frequency", "smoking_status", "moderate_activity")
all_results <- data.frame()
for (i in seq_along(phe_files)) {
    phe_file = file.path("pheno", phe_files[i])
    cat("===== Running:", phe_file, "with", prs_files[i], "=====\n")
    res <- run_logistic(pheno = phe_file, PRS = prs_files[i],
                        covars = covars)
    all_results <- rbind(all_results, res)
}

fwrite(all_results, "NPDs_PRS_logistic_age-sex-centre-BMI-alcohol-etal.csv", sep="\t")









library(ggplot2)
library(dplyr)

df = eur_r2_result
df$Threshold <- factor(df$Threshold, levels = unique(df$Threshold[order(as.numeric(gsub("e-","e-",df$Threshold)))]))

ggplot(df, aes(x = Threshold, y = R2, color = group, group = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(name = expression("McFadden "*R^2)) +
  scale_x_discrete(name = "P-value threshold") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  ggtitle("PRS performance across thresholds")
ggsave("pseudo-r2_eur.png", dpi = 330, width = 8, height = 5)