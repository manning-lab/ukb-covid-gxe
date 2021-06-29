setwd("~/Documents/R/UKB_COVID/FUMA_v2/SNPs/")
COVID_var <- fread("FUMA_snps_v2.txt", stringsAsFactors = F, data.table = F)
setwd("~/Documents/R/UKB_COVID/FUMA_v2/")
HGI <- fread("HGI_UKBB_leave.gz", stringsAsFactors = F, data.table = F)

colnames(HGI)[13] = "rsID"
joint_df <- left_join(COVID_var, HGI, by = "rsID")
joint_df <- joint_df[ ,-c(9:16)]
writexl::write_xlsx(joint_df, path = 'joint_HGI.xlsx')