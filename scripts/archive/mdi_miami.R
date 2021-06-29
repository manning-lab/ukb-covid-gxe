data <- filter(data, data[[pval_col_1]] > 0, data[[pval_col_2]] > 0)
nlps1 <- -log10(mdi[["P_Value_Joint"]])
nlps2 <- -log10(mdi[["P_Value_Interaction"]])


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


top <- cbind(select(mdi, SNP=RSID, CHR, POS), pvalue=nlps1)[keep1, ]
bottom=cbind(select(mdi, SNP=RSID, CHR, POS), pvalue=nlps2)[keep2, ]
topn <- names(top)
bottomn <- names(bottom)
top$Location <- "Top"
bottom$Location <- "Bottom"
d <- rbind(top, bottom)
d$POS <- as.numeric(as.character(d$POS))
d$CHR <- factor(d$CHR, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"))
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
colnames(d_order)[2] <- "Color"

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
            brewer.pal(3, "Dark2")[c(2, 3, 1)])
names(newcols) <-c(levels(factor(lims$Color)), "sex", "cm", "mdi")


d_order <- arrange(d_order, as.integer(Color), desc(pvalue)) %>%
  distinct(SNP, Location, .keep_all=T) 

newcols_bot <-c(rep(x=c("#AAAAAA", "#8A8A8A"), length.out=22, each=1),  # Gray/dark gray for alternating chromosomes
                brewer.pal(3, "Dark2")[c(2, 3, 1)])
names(newcols_bot) <-c(levels(factor(lims$Color)), "sex", "cm", "mdi")


d_order_bot <- arrange(d_order, as.integer(Color), desc(pvalue)) %>%
  distinct(SNP, Location, .keep_all=T) 


#TOP PLOT
p1 <- ggplot() + geom_point(data=d_order[d_order$Location=="Top", ],
                            aes(x=x_coord, y=pvalue, color=factor(Color)),
                            size=0.75, alpha=1) + geom_hline(yintercept=-log10(5e-8), linetype="dashed", color="black") +
  scale_x_continuous(breaks=lims$av[c(1:16, 18, 20, 20, 22)],
                     labels=lims$Color[c(1:16, 18, 20, 20, 22)],
                     expand=c(0,0)) +
  scale_y_continuous(limits=c(2, 8), breaks = seq(2, 8, 2), expand=c(0,0),
                     name=expression(-log[10](italic(p)) - Joint)) +
  scale_colour_manual(name = "Color", values = newcols) +
  scale_fill_manual(name = "Color", values = newcols) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(vjust = -1.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

#BOTTOM PLOT
p2 <- ggplot() + geom_point(data=d_order_bot[d_order_bot$Location== "Bottom",],
                            aes(x=x_coord, y=pvalue, color=factor(Color)), size=0.75, alpha=1) +
  geom_hline(yintercept=-log10(5e-8), linetype="dashed", color="black") +
  scale_x_continuous(breaks=lims$av[c(1:16, 18, 20, 20, 22)],
                     labels=lims$Color[c(1:16, 18, 20, 20, 22)],
                     expand=c(0,0)) +
  scale_y_continuous(limits=c(2, 8), breaks = seq(2, 8, 2), expand=c(0,0),
                     name=expression(-log[10](italic(p)) - Interaction)) +
  scale_colour_manual(name = "Color", values = newcols_bot,
                      # breaks=c(30, 31),
                      # labels=c("Genome-wide significant joint test, no genome-wide significant marginal test",
                      #          "Genome-wide significant interaction test")
  ) +
  scale_fill_manual(name = "Color", values = newcols_bot) +
  #ylab(expression(-log[10](p) - Joint)) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_y_reverse(limits=c(8,2), breaks = seq(8, 2, -2), expand=c(0,0),
                  name=expression(-log[10](italic(p)) - Interaction))


# color_legend <- get_legend(
#   p2 +
#     theme(
#       legend.title=element_blank(),
#       legend.background=element_rect(size=rel(0.5)),
#       legend.key.size=unit(0.1, "lines"),
#       legend.text=element_text(size=10)
#     ) +
#     guides(color=guide_legend(override.aes=list(size=2.5)))
#  )

cowplot::plot_grid(p1 + guides(color=F),
                   p2 + guides(color=F),
                   # color_legend,
                   align='v', nrow=2, rel_heights=c(5, 5))




