setwd("~/Documents/R/UKB_COVID/FUMA_v2")

library(tidyverse)
library(cowplot)
library(data.table)
library(R.utils)
 

sex <- fread("sex_fuma_v2.gz", data.table=F, stringsAsFactors = F)
cm <- fread("cm_fuma_v2.gz", data.table = F, stringsAsFactors = F)
mdi <- fread("mdi_fuma_v2.gz", data.table=F, stringsAsFactors = F)

sex<- transform(sex, CHR = as.character(CHR))
cm<- transform(cm, CHR = as.character(CHR))
mdi<- transform(mdi, CHR = as.character(CHR))


X_sex <- fread("chrX_sex_res_full", data.table=F, stringsAsFactors = F)
X_cm <- fread("chrX_CM_res_full", data.table=F, stringsAsFactors = F)
X_mdi <- fread("chrX_mdi_res_full", data.table=F, stringsAsFactors = F)


xsex<- X_sex[ , c("CHR", "POS", "RSID", "Non_Effect_Allele", "Effect_Allele", "N_Samples", "P_Value_Joint",
                      "P_Value_Interaction", "AF")]

xcm<- X_cm[ , c("CHR", "POS", "RSID", "Non_Effect_Allele", "Effect_Allele", "N_Samples", "P_Value_Joint",
                     "P_Value_Interaction")]

xmdi<- X_mdi[ , c("CHR", "POS", "RSID", "Non_Effect_Allele", "Effect_Allele", "N_Samples", "P_Value_Joint",
                     "P_Value_Interaction")]


### COMBINE THE THREE AND INCLUDE A COLUMN TO KEEPO TRACK OF THE ORIGIN ANALYSIS ###
data <- dplyr::bind_rows(list(sex=sex, cm=cm, mdi=mdi, sex=xsex, cm=xcm, mdi=xmdi), .id="analysis")

### THEN, CONTINUE WITH PROCESSING

# Make a Miami plot -- credit to RaMWAS package (pruning unneeded points)
# and hudson package (basis for Miami plot functionality)
nlps1 <- -log10(data[["P_Value_Joint"]])
nlps2 <- -log10(data[["P_Value_Interaction"]])
# Trim points in crowded regions (credit to RaMWAS package for code snippet)
yfac = as.integer(nlps1 * 100) + 1L
yorder = sort.list(yfac)
levels(yfac) = as.character(seq_len(max(yfac)))
class(yfac) = "factor"
ygroup = split(seq_along(yfac), yfac)
for (i in seq_along(ygroup)) {
  if (length(ygroup[[i]]) > 300) {
    ygroup[[i]] = sample(ygroup[[i]], size=300, replace=FALSE)
  }
}
keep1 = unlist(ygroup, use.names=FALSE)
yfac = as.integer(nlps2 * 100) + 1L
yorder = sort.list(yfac)
levels(yfac) = as.character(seq_len(max(yfac)))
class(yfac) = "factor"
ygroup = split(seq_along(yfac), yfac)
for (i in seq_along(ygroup)) {
  if (length(ygroup[[i]]) > 300) {
    ygroup[[i]] = sample(ygroup[[i]], size=300, replace=FALSE)
  }
}
keep2 = unlist(ygroup, use.names=FALSE)
top <- cbind(dplyr::select(data, SNP=RSID, CHR, POS, analysis), pvalue=nlps1)[keep1, ]
bottom <- cbind(dplyr::select(data, SNP=RSID, CHR, POS, analysis), pvalue=nlps2)[keep2, ]
topn <- names(top)
bottomn <- names(bottom)
top$Location <- "Top"
bottom$Location <- "Bottom"
d <- rbind(top, bottom)
d$POS <- as.numeric(as.character(d$POS))
d$CHR <- factor(d$CHR, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"))
d_order <- d[order(d$CHR, d$POS), ]
d_order$pos_index <- seq.int(nrow(d_order))
chr_lengths <- sapply(1:22, function(chr) max(d[d$CHR == chr, "POS"]))
chr_start_pos <- cumsum(chr_lengths) - chr_lengths
d_order$x_coord <- chr_start_pos[d_order$CHR] + d_order$POS
d_order_sub <- d_order[, c("SNP", "CHR", "POS", "pvalue", "pos_index", "x_coord")]
maxRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.max(x$x_coord),])
minRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.min(x$x_coord),])
milimits <- do.call(rbind, minRows)
malimits <- do.call(rbind, maxRows)
lims <- merge(milimits, malimits, by="CHR")
names(lims) <- c("Color",
                 "snpx", "posx", "px", "posidxx", "xcoordmin",
                 "snpy", "posy", "py", "posidxy", "xcoordmax")
lims$av <- (lims$xcoordmin + lims$xcoordmax)/2
lims <- lims[order(lims$Color),]
# colnames(d_order)[2] <- "Color"
d_order$Color <- ifelse(d_order$pvalue > 6, d_order$analysis, d_order$CHR) # THIS LINE ASSIGNS NEW "COLOR" FOR POINTS ABOVE SOME THRESHOLD
#original; removed
#int_marg_snps <- filter(whr_ss, P_Value_Interaction < 5e-8,P_Value_Marginal < 5e-8)$RSID
#int_nomarg_snps <- filter(whr_ss, P_Value_Interaction < 5e-8,P_Value_Marginal > 5e-8)$RSID
#joint_nomarg_snps <- filter(whr_ss, P_Value_Joint < 5e-8,P_Value_Interaction > 5e-8, P_Value_Marginal > 5e-8)$RSID
# gw_int_snps <- whr_ss$RSID[whr_ss$P_Value_Interaction < 5e-8]
# gw_joint_noMarg_snps <- whr_ss$RSID[whr_ss$P_Value_Joint < 5e-8 &
#                                     whr_ss$P_Value_Main > 5e-8]
#d_order$Color <- case_when(
# d_order$SNP %in% int_marg_snps ~ "30",
#d_order$SNP %in% int_nomarg_snps ~ "31",
#d_order$SNP %in% joint_nomarg_snps ~ "32",
#TRUE ~ as.character(d_order$CHR)
#)
newcols <-c(rep(x=c("#AAAAAA", "#8A8A8A"), length.out=22, each=1),  # Gray/dark gray for alternating chromosomes
            RColorBrewer::brewer.pal(3, "Dark2")[c(2, 3, 1)])
names(newcols) <-c(levels(factor(lims$Color)), "sex", "cm", "mdi")
d_order <- arrange(d_order,Color,
                   desc(pvalue))
                  #%>% distinct(SNP, Location, .keep_all=T)


# newcols_bot <-c(rep(x=c("#AAAAAA", "#8A8A8A"), length.out=22, each=1),  # Gray/dark gray for alternating chromosomes
#                 RColorBrewer::brewer.pal(3, "Dark2")[c(2, 3, 1)])
# names(newcols_bot) <-c(levels(factor(lims$Color)), "sex", "cm", "mdi")
# 
# 
# d_order_bot <- arrange(d_order, as.integer(Color), desc(pvalue)) %>%
#   distinct(SNP, Location, .keep_all=T)

#TOP PLOT
p1 <- ggplot() + geom_point(data=d_order[d_order$Location=="Top", ],
                            aes(x=x_coord, y=pvalue, color=factor(Color)),
                            size=0.75, alpha=1) + geom_hline(yintercept=-log10(5e-8), linetype="dashed", color="black") +
  scale_x_continuous(breaks=lims$av[c(1:16, 18, 20, 20, 22)],
                     labels=lims$Color[c(1:16, 18, 20, 20, 22)],
                     expand=c(0,0)) +
  scale_y_continuous(limits=c(2, 8), breaks = seq(2, 8, 2), expand=c(0,0),
                     name=expression(-log[10](italic(p)) - "Joint Test")) +
  scale_colour_manual(name = "Color", values = newcols, breaks=c("sex", "cm", "mdi"), 
                      labels=c("Sex", "Cardiometabolic", "Multiple Deprivation Index")) +
  #scale_fill_manual(name = "Color", values = newcols) +
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(vjust = -1.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title=element_blank(),
        legend.background=element_rect(color="black", size=rel(1)),
        legend.key.size=unit(0.1, "lines"),
        legend.text=element_text(size=10),
        legend.position = c(.8, .2),
        legend.margin = margin(0,5,5,5))
  
#BOTTOM PLOT
p2 <- ggplot() + geom_point(data=d_order[d_order$Location== "Bottom",],
                            aes(x=x_coord, y=pvalue, color=factor(Color)), size=0.75, alpha=1) +
  geom_hline(yintercept=-log10(5e-8), linetype="dashed", color="black") +
  scale_x_continuous(breaks=lims$av[c(1:16, 18, 20, 20, 22)],
                     labels=lims$Color[c(1:16, 18, 20, 20, 22)],
                     expand=c(0,0)) +
  scale_y_continuous(limits=c(2, 8), breaks = seq(2, 8, 2), expand=c(0,0),
                     name=expression(-log[10](italic(p)) - "Interaction Test")) +
  scale_colour_manual(name = "Color", values = newcols, breaks=c("sex", "cm", "mdi"), 
                      labels=c("Sex", "Cardiometabolic", "Multiple Deprivation Index")
                      # breaks=c(30, 31),
                      # labels=c("Genome-wide significant joint test, no genome-wide significant marginal test",
                      #          "Genome-wide significant interaction test")
  ) +
 #scale_fill_manual(name = "Color", values = newcols_bot) +
  #ylab(expression(-log[10](p) - Joint)) +
 theme_bw()+
   theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_y_reverse(limits=c(8,2), breaks = seq(8, 2, -2), expand=c(0,0),
                  name=expression(-log[10](italic(p)) - "Interaction Test"))
#COLOR LEGEND#
# color_legend <- cowplot::get_legend(
#   p2 + theme(
#       legend.title=element_blank(),
#       legend.background=element_rect(color="black", size=rel(1)),
#       legend.key.size=unit(0.1, "lines"),
#       legend.text=element_text(size=10),
#     )+
#     guides(color=guide_legend(override.aes=list(size=2.5)))
# )
# 
# left_plt <- cowplot::plot_grid(p1 + guides(color=F),
#                    p2 + guides(color=F),
#                    align='v',
#                    #color_legend,
#                    nrow=2, rel_heights=c(5, 5))

cowplot::plot_grid(p1, p2 + guides(color=F), ncol=1, rel_widths = c(1,1))
