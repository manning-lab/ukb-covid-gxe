library(GenomicRanges)
library(tidyverse)

rmarkdown::render("GWIS_pvalue.Rmd", output_file="sexcomp.html",
                   params=list(GEM_file="sex_int_res_full", variant_id_col="RSID",
                  maf_filter= 0.05))
HGI_B2 <- fread("HGI_B2.txt")
HGI_B2 <- filter(HGI_B2, all_inv_var_meta_p<5e-8)
  
sumstats_joint_pruned$start = sumstats_joint_pruned$pos-300000
sumstats_joint_pruned$end = sumstats_joint_pruned$pos+300000

HGI_B2$start = HGI_B2$POS-300000
HGI_B2$end = HGI_B2$POS+300000

gr_sex <- makeGRangesFromDataFrame(sumstats_joint_pruned, seqnames.field = "chr", 
                                   start.field="start", end.field="end")
gr_HGI <- makeGRangesFromDataFrame(HGI_B2, seqnames.field = '#CHR', start.field="start",
                                   end.field="end")
hits = findOverlaps(gr_sex, gr_HGI)
hits


# originally 10 hits and 0 metadata columns...



# Positive Control #

HGI_B2_control <- HGI_B2 %>% add_row(`#CHR`= 1, rsid = "rs12727806", 
                                             start=91700000, end=917800000)
gr_HGI_control <- makeGRangesFromDataFrame(HGI_B2_control, seqnames.field = '#CHR', 
                                           start.field="start", end.field="end")
gr_sex_control <- makeGRangesFromDataFrame(sumstats_joint_pruned, seqnames.field = "chr", 
                                   start.field="start", end.field="end")
hits_control=findOverlaps(gr_sex_control, gr_HGI_control)
hits_control
ranges(gr_sex_control)[queryHits(hits_control)] = ranges(gr_HGI_control)[subjectHits(hits_control)]

# 12 hits and 0 metadata columns

