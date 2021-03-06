library(data.table)
library(tidyverse)
library(R.utils)

sex_int <- fread("leadSNPs_sex_int.txt", stringsAsFactors = F, data.table=F)
sex_joint <- fread("leadSNPs_sexjoint.txt", stringsAsFactors = F, data.table=F)
sex_int <- fread("leadSNPs_sex_int.txt", stringsAsFactors = F, data.table=F)
sex_joint <- fread("leadSNPs_sexjoint.txt", stringsAsFactors = F, data.table=F)
cm_int <- fread("leadSNPs_mdi_int.txt", stringsAsFactors = F, data.table=F)
cm_joint <- fread("leadSNPs_mdi_joint.txt", stringsAsFactors = F, data.table=F)
mdi_int <- fread("leadSNPs_CM_int.txt", stringsAsFactors = F, data.table=F)
mdi_joint <- fread("leadSNPs_cm_joint.txt", stringsAsFactors = F, data.table=F)

cm_int <- cm_int %>% add_column(FUMA = "cm_int")
mdi_int <- mdi_int %>% add_column(FUMA = "mdi_int")
cm_joint <- cm_joint %>% add_column(FUMA = "cm_joint")
mdi_joint <- mdi_joint %>% add_column(FUMA = "mdi_joint")
sex_int <- sex_int %>% add_column(FUMA = "sex_int")
sex_joint <- sex_joint %>% add_column(FUMA = "sex_joint")

fuma_list <- list(cm_int, cm_joint, mdi_int, mdi_joint, sex_int, sex_joint)
fuma_ss <- lapply(fuma_list, function(x) x <- x[,-c(1:3, 8, 9)])
FUMA_SNPs <-do.call("rbind", fuma_ss)
fwrite(FUMA_SNPs,'FUMA_snps.csv',sep='\t')