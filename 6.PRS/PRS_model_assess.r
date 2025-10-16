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



##### PRS 风险分层分析（risk stratification analysis），
library(data.table)
library(dplyr)
library(broom)

risk_stratification  <- function(pheno, PRS){
    # -------------------------
    # 1️⃣ 读取表型数据
    # -------------------------
    phe = fread(pheno)
    phe <- fread(pheno, select = c("FID", "IID", "has_trait"))
    phe$FID = as.character(phe$FID); phe$IID = as.character(phe$IID)
    phe$has_trait <- as.numeric(as.character(phe$has_trait))

    # -------------------------
    # 2️⃣ 加载 PRS 文件并合并
    # -------------------------
    prs <- fread(PRS, select = c("FID", "IID", "SCORE"))
    prs$FID = as.character(prs$FID); prs$IID = as.character(prs$IID)
    phe_prs = merge(phe, prs[, c("FID", "IID", "SCORE")], by=c("FID", "IID"))
    
    phe_prs$scaled = scale(phe_prs$SCORE)

    # 2️按PRS百分位划分
    phe_prs <- phe_prs %>%
    mutate(
        quantile_group = ntile(SCORE, 10),
        top10 = ifelse(quantile_group == 10, 1, 0),
        bottom10 = ifelse(quantile_group == 1, 1, 0)
    )
    
    # 仅保留最高10%与最低10%
    phe_10 <- phe_prs %>% filter(top10 == 1 | bottom10 == 1)
    # 构建模型：高PRS vs 低PRS的风险比
    model <- glm(has_trait ~ top10, data=phe_10, family=binomial)

    # 提取OR与置信区间
    or_table <- exp(cbind(OR = coef(model), confint(model)))
    if ("top10" %in% rownames(or_table)){
        ci_cols <- grep("%", colnames(or_table), value = TRUE)
        or_value <- or_table["top10", "OR"]
        ci_low  <- or_table["top10", ci_cols[1]]
        ci_high <- or_table["top10", ci_cols[2]]
        p_value <- summary(model)$coefficients["top10", 4]
    } else {
        warning("⚠️ top10 variable missing or constant, skipping model")
        or_value <- ci_low <- ci_high <- p_value <- NA
    }

    cat("High vs Low PRS:\n")
    cat("OR =", round(or_value, 3), 
        "95%CI [", round(ci_low, 3), ",", round(ci_high, 3), "]",
        "P =", signif(p_value, 3), "\n\n")

    return(data.frame(
        PHENO = basename(pheno),
        OR = or_value,
        CI_low = ci_low,
        CI_high = ci_high,
        P = p_value
    ))
}

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

all_or <- data.frame()

for (i in seq_along(phe_files)) {
    phe_file = file.path("pheno", phe_files[i])
    cat("===== Running:", phe_file, "=====\n")
    res <- risk_stratification(phe_file, prs_files[i])
    all_or <- rbind(all_or, res)
}

fwrite(all_or, "NPDs_PRS_Top10_vs_Bottom10_OR.csv", sep="\t")