
#R/4.1.0 udunits geos

#library(VennDiagram)
#library(grDevices)
library(RColorBrewer)
library(ggplot2)
library(ggVennDiagram)
library(survival)
library(survminer)
library(ggpubr)
library(dplyr)
library(ggrepel)
library(metafor)

Colors <- read.delim(file = "ColorLibrary", row.names = 1, as.is = T, header = F)

Files <- c("Variant_Matrix.xls.Gene.homo.In_Result.all.corr.xls", "Variant_Matrix.xls.Gene.homo.In_Result.all.CC.corr.xls", "Variant_Matrix.xls.Gene.In_Result.all.corr.xls", "Variant_Matrix.xls.Gene.In_Result.all.CC.corr.xls")

Tabs <- lapply(Files, function(x) {Tab <- read.table(file = x, header = TRUE, row.names = 1, sep = "\t"); names(Tab)[15:16] <- c("variantNum", "GeneLen"); Tab})

names(Tabs) <- c("homoCorr", "homoCCCorr", "hetaCorr", "hetaCCCorr")

TabsNamed <- lapply(names(Tabs), function(x) {names(Tabs[[x]]) <- paste(x, names(Tabs[[x]]), sep = "_"); Tabs[[x]]})


## Meta analysis

Genes <- Reduce(intersect, lapply(Tabs, rownames))

Genes <- Genes[-1]

TabMeta <- matrix(NA, nrow = length(Genes), ncol = 6, dimnames = list(Genes, c("estimate", "se", "ci.lb", "ci.ub", "zval", "pval")))

for(gg in Genes){

	dat <- data.frame(study = names(Tabs),
                  logHR = unlist(lapply(Tabs, function(x) x[gg, "coef"])),
                  selogHR = unlist(lapply(Tabs, function(x) x[gg, "se.coef."])),
                  n = c(385, 224, 385, 224),
                  stringsAsFactors = FALSE)

	## fixed effect

	res1 <- rma(yi = logHR, sei = selogHR, data = dat, method = "FE")

	TabMeta[gg, ] <- c(exp(c(res1$beta, res1$se, res1$ci.lb, res1$ci.ub)), res1$zval, res1$pval)

}

TabMeta <- TabMeta[order(TabMeta[, "pval"], decreasing = F), ]

#lapply(Tabs, function(x) x["MTIF2", ])

#lapply(Tabs, function(x) x["LIMS1", ])

write.table(TabMeta, file = "TabMeta.xls", sep = "\t", quote = F, col.names = NA)

#GeneInfo <- read.delim(file = "../../2.impute/Diff_homo/Variant_Matrix.xls.Gene.Ori.map", as.is = T, header = F)
#
#GeneInfo <- GeneInfo[GeneInfo[[8]]!=".", ]
#
#Tmp <- data.frame(Chr = tapply(GeneInfo[[2]], GeneInfo[[8]], function(x) unique(x)[1]), Start = tapply(GeneInfo[[3]], GeneInfo[[8]], min))
#
##Tmp <- Tmp[which(lapply(Tmp$Chr, length) == 1), ]
#
##Tmp <- Tmp[-1, ]
#
#write.table(Tmp, file = "GeneInfo.xls", quote = F, col.names = NA, sep = "\t")

GeneInfo <- read.delim(file = "GeneInfo.xls", row.names = 1, as.is = T)

ManhattanLis <- list()

for(nn in names(Tabs)){

Tab <- Tabs[[nn]]

Tab <- Tab[rownames(Tab) %in% rownames(GeneInfo), ]

Tab$Chr <- GeneInfo[rownames(Tab), "Chr"]

Tab$Start <- GeneInfo[rownames(Tab), "Start"]


gwasResults <- data.frame(SNP = rownames(Tab), CHR = Tab$Chr, Pos = Tab$Start, P = Tab$Pr...z..)


gwasResults <- gwasResults[with(gwasResults, order(CHR, Pos)), ]

Counts <- table(gwasResults$CHR)

gwasResults$BP <- unlist(lapply(gwasResults$CHR[!duplicated(gwasResults$CHR)], function(x) 1:Counts[x]))

gwasResults$CHR <- as.character(gwasResults$CHR)

gwasResults$CHR <- factor(as.character(gwasResults$CHR), levels = c(as.character(1:22)))

write.table(gwasResults, file = paste0(nn, ".gwasResults.xls"), sep = "\t", quote = F, col.names = NA)

don <- gwasResults %>%

  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%

  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%

  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%

  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

#don$CHR <- factor(as.character(don$CHR), levels = c(as.character(1:22), "X", "Y"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

Inter <- Reduce(intersect, lapply(Tabs[1:4], function(x) {rownames(x)[x$coef > 0 & x$Pr...z.. <= 0.05 & !(is.na(x$Pr...z..))]}))

#donLabel <- don[don$SNP %in% c("LIMS1", "GCC2", "C1orf94", "FHOD3", "PPP2RA", "KRT8P41", "SSB"), ]

Gene23.2 <- unlist(read.table(file = "Gene23", as.is = T))

Gene23 <- c("LIMS1", "GCC2")

#donLabel <- don[don$SNP %in% c("LIMS1", "GCC2", "SIPA1L3", "FHOD3", "C1orf94", "C8orf37-AS1"), ]

donLabel <- don[don$SNP %in% Gene23, ]

donLabel.2 <- don[don$SNP %in% Gene23.2, ]

#stop()

pp1 <- ggplot(don, aes(x=BPcum, y=-log10(P))) +

    # Show all points
    geom_point(aes(color=CHR), alpha=0.8, size=0.6) +

    geom_label_repel(data = donLabel, aes(label = SNP, x=BPcum, y=-log10(P)), min.segment.length = 0) + 

    geom_point(data = donLabel.2, aes(label = SNP, x=BPcum, y=-log10(P)), color = Colors["Red", ], shape = 18, size = 3) +

    scale_color_manual(values = rep(Colors[c("Blue2", "Orange"), ], 22 )) +

    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(NA, 4)) +     # remove space between plot area and x axis

    # Custom the theme:

    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = Colors["Red", ], size = 0.3) +

    xlab("Chromosome") +

    theme_bw() +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      plot.margin = margin(10, 10, 20, 10, "pt")
    )

ManhattanLis[[nn]] <- pp1

}

ggsave("ManhattanPlots.pdf", ggarrange(ManhattanLis[["homoCorr"]], ManhattanLis[["homoCCCorr"]], ManhattanLis[["hetaCorr"]], ManhattanLis[["hetaCCCorr"]], nrow = 4, ncol = 1), width = 14, height = 12)


pp <- ggVennDiagram(lapply(Tabs[1:4], function(x) {rownames(x)[x$coef > 0 & x$Pr...z.. <= 0.01 & !(is.na(x$Pr...z..))]}), label_alpha = 0, category.names = c("Double","Double \nE->E","Any", "Any \nE->E"), label = "count", label_size = 4, set_size = 4, color = "black") + theme(legend.position = "null", plot.margin = margin(10, 10, 10, 10, "pt")) + scale_fill_gradient(low=brewer.pal(9, "YlOrRd")[2],high = brewer.pal(9, "YlOrRd")[6])


pp05 <- ggVennDiagram(lapply(Tabs[1:4], function(x) {rownames(x)[x$coef > 0 & x$Pr...z.. <= 0.05 & !(is.na(x$Pr...z..))]}), label_alpha = 0, category.names = c("Double","Double \nE->E","Any", "Any \nE->E"), label = "count", label_size = 4, set_size = 4, color = "black") + theme(legend.position = "null", plot.margin = margin(0, 0, 0, 0)) + scale_fill_gradient(low=brewer.pal(9, "YlOrRd")[2],high = brewer.pal(9, "YlOrRd")[6])

ggsave("Venn_005.pdf", pp05)

#pdf(file="venn_0.01_Corr.pdf")
#
#        temp <- venn.diagram(lapply(Tabs[1:4], function(x) {rownames(x)[x$coef > 0 & x$Pr...z.. <= 0.01 & !(is.na(x$Pr...z..))]}), fill = brewer.pal(n = 11, name = "Set3")[c(1,2,4,5)], alpha = 1, filename = NULL)
#
#        grid.draw(temp)
#
#dev.off()

Inter <- Reduce(intersect, lapply(Tabs[1:4], function(x) {rownames(x)[x$coef > 0 & x$Pr...z.. <= 0.05 & !(is.na(x$Pr...z..))]}))

write.table(Inter, file = "InterCorr05.list", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(Reduce(cbind, lapply(TabsNamed[1:4], function(x) x[Inter, ])), file = "InterCorr0.05.xls", quote = FALSE, sep = "\t", col.names = NA)

Tmp <- Reduce(cbind, lapply(TabsNamed[1:4], function(x) x[Inter, c(2,5,6,15,16)])) 

Tmp$Chr <- GeneInfo[rownames(Tmp), "Chr"]

write.table(Tmp, file = "InterCorr0.05.sim.xls", quote = FALSE, sep = "\t", col.names = NA)


GOCAR <- read.delim(file = "GOCAR.xls", row.names = 1, as.is = T)

GOCAR$graft_loss_time_update <- GOCAR$graft_loss_time_update/30

Heta <- read.delim(file = "Variant_Matrix.xls.Gene.In", row.names = 1, as.is = T, check.names = F)

Homo <- read.delim(file = "Variant_Matrix.xls.Gene.homo.In", row.names = 1, as.is = T, check.names = F)

GOCAR$LIMS1_Heta_freq <- unlist(Heta["LIMS1", rownames(GOCAR)])

GOCAR$LIMS1_Heta <- GOCAR$LIMS1_Heta_freq/IQR(GOCAR$LIMS1_Heta_freq)

GOCAR$LIMS1_Homo_freq <- unlist(Homo["LIMS1", rownames(GOCAR)])

GOCAR$LIMS1_Homo <- GOCAR$LIMS1_Homo_freq

GOCAR$R2R <- paste(GOCAR$Genetic_Donor_Race, GOCAR$Genetic_Rec_Race, sep = " > ")

GOCARCC <- GOCAR[GOCAR$R2R == "Caucasian > Caucasian", ]


#stop()

GOCAR$LIMS1_HetaG <- cut(GOCAR$LIMS1_Heta, breaks = c(min(GOCAR$LIMS1_Heta)- 0.1, mean(GOCAR$LIMS1_Heta), max(GOCAR$LIMS1_Heta)))

levels(GOCAR$LIMS1_HetaG) <- c("Low", "High")

fit <- survfit(Surv(graft_loss_time_update, death_censored_graft_loss_event) ~ LIMS1_HetaG, data = GOCAR)

pp2 <- ggsurvplot(fit, palette = Colors[c("Blue", "Red"), ], ggtheme = theme_classic(), conf.int = F, pval = TRUE, pval.size = 4, legend.labs=c("Lower \nquantile", "Upper \nquantile"), font.x = 10, font.y = 10, font.legend = 10, font.tickslab = 10, legend.title = "LIMS1 \nmismatch (any)", xlab = "Follow-up time (months)", risk.table = T, risk.table.col = "strata", font.title = c(10, "bold", "darkblue"), risk.table.fontsize = 4, risk.table.y.text = F)

pp2$plot <- pp2$plot + theme(plot.margin = margin(0, 25, 0, 18, "pt"))

pp2$table <- pp2$table +
  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")


GOCAR$LIMS1_HomoG <- cut(GOCAR$LIMS1_Homo, breaks = c(min(GOCAR$LIMS1_Homo)- 0.1, mean(GOCAR$LIMS1_Homo), max(GOCAR$LIMS1_Homo)))

levels(GOCAR$LIMS1_HomoG) <- c("Low", "High")

fit <- survfit(Surv(graft_loss_time_update, death_censored_graft_loss_event) ~ LIMS1_HomoG, data = GOCAR)

pp3 <- ggsurvplot(fit, palette = Colors[c("Blue", "Red"), ], ggtheme = theme_classic(), conf.int = F, pval = TRUE, pval.size = 4, legend.labs=c("Lower \nquantile", "Upper \nquantile"), font.x = 10, font.y = 10, font.legend = 10, font.tickslab = 10, legend.title = "LIMS1 \nmismatch (double)", xlab = "Follow-up time (months)", risk.table = T, risk.table.col = "strata", risk.table.y.text = F, font.title = c(10, "plain", "darkblue"), risk.table.fontsize = 4)

pp3$plot <- pp3$plot + theme(plot.margin = margin(0, 25, 0, 18, "pt"))

pp3$table <- pp3$table +
	  theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 20, 0, 8, "pt")) + ylab("") + xlab("")

#ggsave("FigLIMS1.pdf", ggarrange(ManhattanLis[["hetaCorr"]], ggarrange(pp, ggarrange(pp2$plot, pp2$table, nrow = 2, heights = c(2, 1)), ggarrange(pp3$plot, pp3$table, nrow = 2, heights = c(2, 1)), labels = c("B", "C", "D"), ncol = 3, nrow = 1, widths = c(3.5, 2, 2)), nrow = 2, labels = "A"), width = 15, height = 10)

ggsave("FigLIMS1.pdf", ggarrange(ManhattanLis[["hetaCorr"]], ggarrange(pp, ggarrange(pp2$plot, pp2$table, nrow = 2, heights = c(2, 1)), ggarrange(pp3$plot, pp3$table, nrow = 2, heights = c(2, 1)), ncol = 3, nrow = 1, widths = c(3.5, 2, 2)), nrow = 2), width = 15, height = 10)

#ggsave("Test.pdf", pp2$plot)


GOCARCC$LIMS1_HetaG <- cut(GOCARCC$LIMS1_Heta, breaks = c(min(GOCARCC$LIMS1_Heta)- 0.1, mean(GOCARCC$LIMS1_Heta), max(GOCARCC$LIMS1_Heta)))

levels(GOCARCC$LIMS1_HetaG) <- c("Low", "High")

fit <- survfit(Surv(graft_loss_time_update, death_censored_graft_loss_event) ~ LIMS1_HetaG, data = GOCARCC)

pp2C <- ggsurvplot(fit, palette = Colors[c("Blue", "Red"), ], ggtheme = theme_classic(), conf.int = F, pval = TRUE, pval.size = 4, legend.labs=c("Lower \nquantile", "Upper \nquantile"), font.x = 10, font.y = 10, font.legend = 10, font.tickslab = 10, legend.title = "LIMS1 \nmismatch (any)", xlab = "Follow-up time (months)", risk.table = T, font.title = c(10, "bold", "darkblue"), risk.table.col = "strata", risk.table.y.text = F, risk.table.fontsize = 4)

pp2C$plot <- pp2C$plot + theme(plot.margin = margin(0, 5, 0, 18))
    
pp2C$table <- pp2C$table +
 	 theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none") + ylab("") + xlab("")

GOCARCC$LIMS1_HomoG <- cut(GOCARCC$LIMS1_Homo, breaks = c(min(GOCARCC$LIMS1_Homo)- 0.1, mean(GOCARCC$LIMS1_Homo), max(GOCARCC$LIMS1_Homo)))

levels(GOCARCC$LIMS1_HomoG) <- c("Low", "High")

fit <- survfit(Surv(graft_loss_time_update, death_censored_graft_loss_event) ~ LIMS1_HomoG, data = GOCARCC)

pp3C <- ggsurvplot(fit, palette = Colors[c("Blue", "Red"), ], ggtheme = theme_classic(), conf.int = F, pval = TRUE, pval.size = 4, legend.labs=c("Lower \nquantile", "Upper \nquantile"), font.x = 10, font.y = 10, font.legend = 10, font.tickslab = 10, legend.title = "LIMS1 \nmismatch (double)", xlab = "Follow-up time (months)", risk.table = T, font.title = c(10, "plain", "darkblue"), risk.table.fontsize = 4, risk.table.col = "strata", risk.table.y.text = F)

pp3C$plot <- pp3C$plot + theme(plot.margin = margin(0, 5, 0, 18))

pp3C$table <- pp3C$table +
          theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none") + ylab("") + xlab("")


ggsave("FigLIMS1_S.pdf", ggarrange(pp05, ggarrange(pp2C$plot, pp2C$table, nrow = 2, heights = c(2, 1)), ggarrange(pp3C$plot, pp3C$table, nrow = 2, heights = c(2, 1)), labels = c("A", "B", "C"), ncol = 3, nrow = 1, widths = c(3.5, 2, 2)), width = 15, height = 5)

