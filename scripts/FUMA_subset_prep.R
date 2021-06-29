library(data.table)
library(tidyverse)
library(R.utils)

sex_full[which(sex_full$RSID == "rs182113773"),]
mdi_full[which(mdi_full$RSID == "rs576999056"),]

var <- HGI[which(HGI$rsid == "rs148793499"), ]
rs6747163 <- BMI_T2D[which(BMI_T2D$RSID =="rs6747163"),]
rs2268616 <- sex[which(sex$RSID == "rs2268616"), ]
var <- HGI[which(HGI$rsid == "rs188530487"), ]
var2 <- sex[which(sex$POS == 218260234), ]
rs11115199 <- cm[which(cm$RSID== "rs11115199"), ]
rs148793499 <- cm[which(cm$RSID== "rs148793499"), ]

rs12461506 <- sex[which(sex$RSID =="rs12461506"),]
rs114807731 <- sex[which(sex$RSID =="rs114807731"),]

rs148817892 <- mdi[which(mdi$RSID =="rs148817892"),]
rs117488928 <- mdi[which(mdi$RSID =="rs117488928"),]
rs192911167 <- mdi[which(mdi$RSID =="rs192911167"),]
rs113541905 <- mdi[which(mdi$RSID =="rs113541905"),]
rs11945368 <- mdi[which(mdi$RSID =="rs11945368"),]

setwd("~/Documents/R/UKB_COVID/FUMA_v2/")
sex <- fread("sex_full.txt.gz", stringsAsFactors = F, data.table=F)
sex_sb <-sex[c("CHR", "POS", "RSID", "Non_Effect_Allele", "Effect_Allele", 
               "N_Samples", "P_Value_Joint", "P_Value_Interaction", "AF")]
fwrite(sex_sb,'sex_fuma_v2',sep='\t')
gzip('sex_fuma_v2',destname='sex_fuma_v2.gz')
       
                             
mdi <- fread("mdi_full.txt.gz", stringsAsFactors = F, data.table=F)
mdi_sb <-mdi[c("CHR", "POS", "RSID", "Non_Effect_Allele", "Effect_Allele", 
                         "N_Samples", "P_Value_Joint", "P_Value_Interaction", "AF")]
fwrite(mdi_sb,'mdi_fuma',sep='\t')
gzip('mdi_fuma',destname='mdi_fuma_v2.gz')

cm <- fread("joint_bmi_t2d_full.txt.gz", stringsAsFactors = F, data.table=F)
cm_sb <-cm[c("CHR", "POS", "RSID", "Non_Effect_Allele", "Effect_Allele", 
             "N_Samples", "P_Value_Joint", "P_Value_Interaction", "AF")]
fwrite(cm_sb,'cm_fuma',sep='\t')
gzip('cm_fuma',destname='cm_fuma_v2.gz')



sex_sb <- arrange(sex_sb, P_Value_Joint)
sex_sb <- head(sex_sb, 20)
write_xlsx(sex_sb, path = 'sex_full.xlsx')
