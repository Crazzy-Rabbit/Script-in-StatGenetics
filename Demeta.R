#=========================================================================#
#                      de meta for Meta and AJHG
#=========================================================================#
rescale_b_se <- function(gwas, case, control){
  names(gwas) = c("SNP","A1","A2","freq","b","se","P","N")
  # n           = gwas$N
  # i           = dnorm(qnorm(1-K))/K
  # n_eff       = i^2*v*(1-v)/(1-K)^2 * n # Yang et al. Gen Epi 2010
  n_eff       = as.integer(4*case*control/(case + control))
  p           = gwas$freq
  z           = gwas$b / gwas$se
  se          = 1/sqrt(2*p*(1-p)*(n_eff+z^2))
  b           = z*se
  gwas$b      = b
  gwas$se     = se
  gwas$N      = n_eff
  
  return(gwas)
}

de_meta <- function(meta, sub_c){
  if (!requireNamespace("logger", quietly=TRUE)){ install.packages("logger") }
  lapply(c("dplyr","data.table","logger"), function(pkg) suppressWarnings(library(pkg,character.only=TRUE)))
  
  log_info("Start to check the input meta...")
  log_info("Start to check SNP effect size")
  if(any(meta$beta_meta == 0)){
    Delmeta = meta[meta$beta_meta == 0, ]
    log_info(paste("The meta file contained", nrow(Delmeta), "SNPs that effect size = 0"))
    meta = meta[!(SNP %in% Delmeta$SNP)]
    log_info("We will delete these SNPs")
  }else{ log_info("No SNPs effect size = 0") }
  
  log_info("Start to check P value")
  if(any(meta$p > 1 | meta$p < 0)){
    Delmeta = meta[meta$p > 1 | meta$p < 0, ]
    log_info(paste("The meta file contained", nrow(Delmeta), "SNPs that Pvalue > 1 or < 0"))
    meta = meta[!(SNP %in% Delmeta$SNP)]
    log_info("We will delete these SNPs")
  }else{ log_info("No SNPs Pvalue > 1 or < 0") }
  
  log_info("Start to check the input sub_c...")
  log_info("Start to check SNP effect size")
  if(any(sub_c$beta_c == 0)){
    Delmeta = sub_c[sub_c$beta_c == 0, ]
    log_info(paste("The sub_c file contained", nrow(Delmeta), "SNPs that effect size = 0"))
    sub_c = sub_c[!(SNP %in% Delmeta$SNP)]
    log_info("We will delete these SNPs")
  }else{ log_info("No SNPs effect size = 0") }
  
  log_info("Start to merge the meta and sub_c by SNP")
  combined_data <- merge(meta,sub_c,by="SNP")
  
  log_info("Start to keep meta only SNPs")
  meta_only <- meta[meta$SNP %in% setdiff(meta$SNP,sub_c$SNP), ]
  meta_only [, n_meta := n_meta - unique(sub_c$N)]
  setDT(meta_only)
  setnames(meta_only, c("SNP","A1","A2","freq","b","se","p","N"))
  log_info(paste("There were", nrow(meta_only), "SNPs only occured in meta"))
  
  log_info("Start to fix the effect allele in meta and sub_c")
  combined_data[, beta_c := ifelse((A1_m==A1_c & A2_m==A2_c),beta,ifelse(A1_m==A2_c & A2_m==A1_c,-beta,NA_real_))]
  if(any(is.na(combined_data$beta_c))) {
    Deldt = subset(combined_data, is.na(combined_data$beta_c))
    log_info(paste("There were", nrow(Deldt), "SNPs have unmatched alleles, will delete these SNPs!"))
    log_info(paste("In these SNPs, the lowest P value is", min(Deldt$p)))
    combined_data = combined_data[!(SNP %in% Deldt$SNP)]
  }
  
  log_info("Start to demeta SNP effect size and se")
  se_meta_mean <- mean(combined_data$se_meta, na.rm=TRUE)
  combined_data[, se_demeta := fifelse((1/se_meta^2 - 1/se_c^2) > 0, sqrt(1/(1/se_meta^2 - 1/se_c^2)),se_meta_mean)]
  combined_data[, beta_demeta := beta_meta-se_demeta^2 * (beta_c-beta_meta)/se_c^2]
  
  
  if (any(combined_data$beta_demeta == 0)) {
    log_info("There may be generate zero beta after demeta!")
    Deldt = combined_data[combined_data$beta_demeta == 0, ]
    log_info(paste("There were", nrow(Deldt), "SNPs has been deleted casued by the demeta value is zero!"))
    combined_data = combined_data[!(SNP %in% Deldt$SNP), ]
  }
  
  log_info("Start to create effect N of demeta")  
  combined_data[, N_eff := n_meta - N]
  
  common_SNP <- combined_data[, c("SNP","A1_m","A2_m","freq","beta_demeta","se_demeta","p","N_eff")]
  setDT(common_SNP)
  setnames(common_SNP, c("SNP","A1","A2","freq","b","se","p","N"))
  
  log_info("Start to combine meta only SNPs and demeta SNPs to one")  
  outfile <- rbind(common_SNP,meta_only)
  
  log_info("Start to calculate P value based on beta and se")
  outfile$p <- 2 * (1 - pnorm(abs(outfile$b / outfile$se)))
  
  log_info("Demeta analysis is completed!") 
  return(outfile)
}

library(data.table)
setwd("E:/Shi/ALL_of_my_Job/24-28川大华西/2_project_hearing loss/process/PRS/PRS_UKB/demeta")
hm3 = fread("listHM3.txt", col.names="SNP")

meta = fread("ARHL_MVP_AJHG_BBJ_reformatMETAL.txt")
meta_hm3 = merge(meta, hm3, by="SNP")
names(meta_hm3) = c("SNP","A1_m","A2_m","freq","beta_meta","se_meta","p","n_meta")

ajhg = fread("AJHG_EUR_reformat.gctb.txt")
ajhg_hm3 = merge(ajhg, hm3, by="SNP")
names(ajhg_hm3) = c("SNP","A1_m","A2_m","freq","beta_meta","se_meta","p","n_meta")

## rescale beta se for sub cohort data
sub_c = fread("2247_1.v1.1.fastGWA", select=c("SNP","A1","A2","AF1","BETA","SE","P","N"))
sub_c_hm3 = merge(sub_c, hm3, by="SNP")

re_sub_c_hm3 = rescale_b_se(sub_c_hm3, 114318, 323449)
names(re_sub_c_hm3) = c("SNP","A1_c","A2_c","freq_c","beta","se_c","p_c","N")

# run demeta for meta and ajhg
demeta_meta_hm3 = de_meta(meta_hm3, re_sub_c_hm3)
demeta_ajhg_hm3 = de_meta(ajhg_hm3, re_sub_c_hm3)

ggplot(demeta_meta_hm3, aes(x = N)) +
  geom_histogram(fill = "#69b3a2", color = "#e9ecef", bins = 100, alpha = 0.8) +
  scale_x_continuous(labels = scales::comma) +
  labs(title="Distribution of Sample Size (N)", x="Sample Size", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggplot(demeta_ajhg_hm3, aes(x = N)) +
  geom_histogram(fill = "#69b3a2", color = "#e9ecef", bins = 100, alpha = 0.8) +
  scale_x_continuous(labels = scales::comma) +
  labs(title = "Distribution of Sample Size (N)",x = "Sample Size", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

out_meta_demeta_hm3 = demeta_meta_hm3[which(demeta_meta_hm3$N > 100000),]
out_ajhg_demeta_hm3 = demeta_ajhg_hm3[which(demeta_ajhg_hm3$N > 100000),]
fwrite(out_meta_demeta_hm3, file="Demeta_raw_hm3_AJHG_rescale_keep_metaonly_SNP_07May2025.txt", sep="\t")
fwrite(out_ajhg_demeta_hm3, file="Demeta_raw_hm3_Meta_rescale_keep_metaonly_SNP_07May2025.txt", sep="\t")