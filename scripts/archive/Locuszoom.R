library(data.table)
library(tidyverse)
library(R.utils)

setwd("~/Documents/R/UKB_COVID/FUMA_v2/")
sex <- fread("sex_full.txt.gz", stringsAsFactors = F, data.table=F) #GEM File
BMI_T2D <- fread("joint_bmi_t2d_full.txt.gz", stringsAsFactors = F, data.table=F)
mdi <- fread("mdi_full.txt.gz", stringsAsFactors = F, data.table=F)
HGI <- fread("HGI_UKBB_leave.gz", stringsAsFactors = F, data.table = F)

#rs11115199 = 12:82510665, cm_int ---- chr 14 instead of 12 for GEM, HGI good
#Rs148793499 -  18:58314588, cm_joint ---- chr 2 instead of 18 for GEM
#Rs2268616 - 14:75419444, sexjoint --- good 
#2:21820234, sex --- good
var1 <- BMI_T2D[which(BMI_T2D$RSID == "rs11115199"), ]
var2 <- BMI_T2D[which(BMI_T2D$RSID == "rs148793499"), ]


var3 <- sex[which(sex$RSID == "rs2268616"), ]
var4 <- sex[which(sex$RSID == "2:218260234_AC_A"), ]

gw_loci <- rbind(var1, var2)
gw_loci_sex <- rbind(var3, var4)

for (i in 1:nrow(gw_loci)) {
  gw_chr <- gw_loci$CHR[i]
  gw_pos <- gw_loci$POS[i]
  lz_track <- BMI_T2D %>%
    mutate(CHR=as.integer(CHR)) %>%
    filter(CHR == gw_chr,
           POS > gw_pos - 500000, 
           POS < gw_pos + 500000) %>%
    select(CHR, start=POS, end=POS, rsID=RSID, SNP=SNPID, REF=Non_Effect_Allele, ALT=Effect_Allele, 
           AF, Beta=Beta_Interaction_1, Pmain=P_Value_Marginal, 
           Pint=P_Value_Interaction, Pjoint=P_Value_Joint)
  int_fname <- paste0("../../UKB_COVID/Locuszoom/", 
                      gsub(":", "_", gw_loci$SNPID[i]), "_track.tsv")
  write_tsv(lz_track, int_fname)
  system(paste0("bgzip -f ", int_fname))
  system(paste0("tabix -f -p bed -S 1 ", int_fname, ".gz"))
}

#HGI 
var1_hgi <- HGI[which(HGI$rsid == "rs11115199"), ]
var2_hgi <- HGI[which(HGI$rsid == "rs9946577"), ]
#proxy for rs148793499 (18:58314588) to create 50 kb range
# rs9946577            (18:58314520)
var3_hgi <- HGI[which(HGI$rsid == "rs2268616"), ]
var4_hgi <- HGI[which(HGI$rsid == "rs2061808"), ]
#proxy for 2:218260234 to create 50 kb range
#rs2061808 2:218259884

gw_loci <- rbind(var1_hgi, var2_hgi, var3_hgi, var4_hgi)

for (i in 1:nrow(gw_loci)) {
  gw_chr <- gw_loci$`#CHR`[i]
  gw_pos <- gw_loci$POS[i]
  lz_track <- HGI %>%
    mutate(CHR=as.integer(`#CHR`)) %>%
    filter(CHR == gw_chr,
           POS > gw_pos - 500000, 
           POS < gw_pos + 500000) %>%
    select(CHR=`#CHR`, start=POS, end=POS, rsID=rsid, SNP, REF, ALT, 
           AF=all_meta_AF, Beta=all_inv_var_meta_beta, Pmain=all_inv_var_meta_p)
  int_fname <- paste0("../../UKB_COVID/Locuszoom/", 
                      gsub(":", "_", gw_loci$SNP[i]), "_track.tsv")
  write_tsv(lz_track, int_fname)
  system(paste0("bgzip -f ", int_fname))
  system(paste0("tabix -f -p bed -S 1 ", int_fname, ".gz"))
}




#####################
    
cm <- fread("sig_cm_var.csv", stringsAsFactors = F, data.table = F)
mdi <- fread("sig_mdi_var.csv", stringsAsFactors = F, data.table = F)
sex <- fread("sig_sex_var.csv", stringsAsFactors = F, data.table = F)


chr18_hgi <- HGI[which(HGI$`#CHR` == "2"), ]
pos_hgi <- chr18_hgi[which(chr18_hgi$POS < 218260234 + 500), ]
pos_hgi <- pos_hgi[which(pos_hgi$POS > 218260234 - 500), ]
