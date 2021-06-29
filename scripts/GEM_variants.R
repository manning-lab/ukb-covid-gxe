library(data.table)
library(tidyverse)
library(R.utils)

sig_var <- readxl::read_excel("UKBB_sig_var.xlsx") #FUMA GWAS sig variants
colnames(sig_var)[1] = "RSID"
sex <- fread("sex_full.txt.gz", stringsAsFactors = F, data.table=F) #GEM File

sex_joined <- left_join(sig_var, sex, by = "RSID")
sex_joined <- sex_joined[,-c(6:8)]
colnames(sex_joined)[8] = "GEM_n"
colnames(sex_joined)[14] = "GEM_p_marginal"
colnames(sex_joined)[15] = "GEM_p_interaction"
colnames(sex_joined)[16] = "GEM_p_joint"
sex_joined <- sex_joined[-c(21),]
fwrite(sex_joined,'sig_sex_var.csv',sep='\t')

BMI_T2D <- fread("joint_bmi_t2d_full.txt.gz", stringsAsFactors = F, data.table=F)
BMI_joined <- left_join(sig_var, BMI_T2D, by = "RSID")
BMI_joined <- BMI_joined[,-c(6:8)]
colnames(BMI_joined)[8] = "GEM_n"
colnames(BMI_joined)[18] = "GEM_p_marginal"
colnames(BMI_joined)[19] = "GEM_p_interaction"
colnames(BMI_joined)[20] = "GEM_p_joint"
BMI_joined <- BMI_joined[-c(21),]
fwrite(BMI_joined,'sig_cm_var.csv',sep='\t')

mdi <- fread("mdi_full.txt.gz", stringsAsFactors = F, data.table=F)
mdi_joined <- left_join(sig_var, mdi, by = "RSID")
mdi_joined <- mdi_joined[,-c(6:8)]
colnames(mdi_joined)[8] = "GEM_n"
colnames(mdi_joined)[14] = "GEM_p_marginal"
colnames(mdi_joined)[15] = "GEM_p_interaction"
colnames(mdi_joined)[16] = "GEM_p_joint"
mdi_joined <- mdi_joined[-c(21),]
fwrite(mdi_joined,'sig_mdi_var.csv',sep='\t')


HGI <- fread("HGI_UKBB_leave.gz", stringsAsFactors = F, data.table = F)
proxies <- fread("COVID_UKBB_rsIDs_proxies.txt", stringsAsFactors = F, data.table = F)
var <- HGI[which(HGI$rsid == "rs73162794"), ]
