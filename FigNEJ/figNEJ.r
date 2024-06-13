
#R/4.1.0 udunits geos

library(RColorBrewer)
library(ggplot2)
library(ggVennDiagram)
library(survival)
library(survminer)
library(ggpubr)

GOCAR <- read.delim(file = "GOCAR.xls", row.names = 1, as.is = T)

GOCAR$graft_loss_time_update <- GOCAR$graft_loss_time_update/30

GOCAR$NEJGenoGroup <- factor(GOCAR$NEJGenoGroup, levels = c("Other", "NoRisk", "Risk"))


Colors <- read.delim(file = "ColorLibrary", row.names = 1, as.is = T, header = F)

fit <- survfit(Surv(graft_loss_time_update, death_censored_graft_loss_event) ~ NEJRiskHeta, data = GOCAR)

pp2 <- ggsurvplot(fit, palette = Colors[c("Blue", "Red"), ], ggtheme = theme_classic(), conf.int = F, pval = TRUE, pval.size = 4, legend.labs=c("Non-risk", "Risk"), font.x = 10, font.y = 10, font.legend = 10, font.tickslab = 10, legend.title = "rs893403 risk\nmismatch", xlab = "Follow-up time (months)", risk.table = T, font.title = c(10, "bold", "darkblue"), risk.table.col = "strata", risk.table.y.text = F, risk.table.fontsize = 4)

pp2$plot <- pp2$plot + theme(plot.margin = margin(0, 5, 0, 18))

pp2$table <- pp2$table +
          theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none") + ylab("") + xlab("")


fit <- survfit(Surv(graft_loss_time_update, death_censored_graft_loss_event) ~ NEJGenoGroup, data = GOCAR)

pp3 <- ggsurvplot(fit, palette = Colors[c("Blue", "Yellow", "Red"), ], ggtheme = theme_classic(), conf.int = F, pval = TRUE, pval.size = 4, legend.labs=c("Other", "Non-risk", "Risk"), font.x = 10, font.y = 10, font.legend = 10, font.tickslab = 10, legend.title = "rs893403 risk\nmismatch", xlab = "Follow-up time (months)", risk.table = T, font.title = c(10, "bold", "darkblue"), risk.table.col = "strata", risk.table.y.text = F, risk.table.fontsize = 4)

pp3$plot <- pp3$plot + theme(plot.margin = margin(0, 5, 0, 18))

pp3$table <- pp3$table +
          theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none") + ylab("") + xlab("")


fit <- survfit(Surv(graft_loss_time_update, death_censored_graft_loss_event) ~ NEJGenoGroup, data = GOCAR[GOCAR$R2R == "Caucasian > Caucasian", ])

pp3.cc <- ggsurvplot(fit, palette = Colors[c("Blue", "Yellow", "Red"), ], ggtheme = theme_classic(), conf.int = F, pval = TRUE, pval.size = 4, legend.labs=c("Other", "Non-risk", "Risk"), font.x = 10, font.y = 10, font.legend = 10, font.tickslab = 10, legend.title = "rs893403 risk\nmismatch", xlab = "Follow-up time (months)", risk.table = T, font.title = c(10, "bold", "darkblue"), risk.table.col = "strata", risk.table.y.text = F, risk.table.fontsize = 4)

pp3.cc$plot <- pp3.cc$plot + theme(plot.margin = margin(0, 5, 0, 18))

pp3.cc$table <- pp3.cc$table +
          theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none") + ylab("") + xlab("")


fit <- survfit(Surv(graft_loss_time_update, death_censored_graft_loss_event) ~ NEJGenoGroup, data = GOCAR[GOCAR$R2R != "Caucasian > Caucasian", ])

pp3.nocc <- ggsurvplot(fit, palette = Colors[c("Blue", "Yellow", "Red"), ], ggtheme = theme_classic(), conf.int = F, pval = TRUE, pval.size = 4, legend.labs=c("Other", "Non-risk", "Risk"), font.x = 10, font.y = 10, font.legend = 10, font.tickslab = 10, legend.title = "rs893403 risk\nmismatch", xlab = "Follow-up time (months)", risk.table = T, font.title = c(10, "bold", "darkblue"), risk.table.col = "strata", risk.table.y.text = F, risk.table.fontsize = 4)

pp3.nocc$plot <- pp3.nocc$plot + theme(plot.margin = margin(0, 5, 0, 18))

pp3.nocc$table <- pp3.nocc$table +
          theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none") + ylab("") + xlab("")


fit <- survfit(Surv(graft_loss_time_update, death_censored_graft_loss_event) ~ NEJGenoGroup, data = GOCAR)

pp3 <- ggsurvplot(fit, palette = Colors[c("Blue", "Yellow", "Red"), ], ggtheme = theme_classic(), conf.int = F, pval = TRUE, pval.size = 4, legend.labs=c("Other", "Non-risk", "Risk"), font.x = 10, font.y = 10, font.legend = 10, font.tickslab = 10, legend.title = "rs893403 risk\nmismatch", xlab = "Follow-up time (months)", risk.table = T, font.title = c(10, "bold", "darkblue"), risk.table.col = "strata", risk.table.y.text = F, risk.table.fontsize = 4)

pp3$plot <- pp3$plot + theme(plot.margin = margin(0, 5, 0, 18))

pp3$table <- pp3$table +
          theme(plot.title = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = "none") + ylab("") + xlab("")

ggsave("FigNEJDirection.pdf", ggarrange(ggarrange(pp3$plot, pp3$table, nrow = 2, heights = c(2, 1)), ggarrange(pp3.cc$plot, pp3.cc$table, nrow = 2, heights = c(2, 1)), ggarrange(pp3.nocc$plot, pp3.nocc$table, nrow = 2, heights = c(2, 1)), ncol = 3), height = 5, width = 12)


        #Clinical$GenoGroupNEJHaplo <- factor(Clinical$GenoGroupNEJHaplo, levels = c("Other", "Risk", "Min", "Maj"))

        fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", "NEJGenoGroup", " + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum")), data = GOCAR)

        Result <- cbind(paste0("(", round(summary(fit)$conf.int[, 3, drop = FALSE], 2), ", ", round(summary(fit)$conf.int[, 4, drop = FALSE], 2), ")"), round(summary(fit)$coefficients[, c(2, 5)], 2), summary(fit)$n)

        fit.cc <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", "NEJGenoGroup", " + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum")), data = GOCAR[GOCAR$R2R == "Caucasian > Caucasian", ])

        Result.cc <- cbind(paste0("(", round(summary(fit.cc)$conf.int[, 3, drop = FALSE], 2), ", ", round(summary(fit.cc)$conf.int[, 4, drop = FALSE], 2), ")"), round(summary(fit.cc)$coefficients[, c(2, 5)], 2), summary(fit.cc)$n)

        fit.nocc <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", "NEJGenoGroup", " + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum")), data = GOCAR[GOCAR$R2R != "Caucasian > Caucasian", ])

        Result.nocc <- cbind(paste0("(", round(summary(fit.nocc)$conf.int[, 3, drop = FALSE], 2), ", ", round(summary(fit.nocc)$conf.int[, 4, drop = FALSE], 2), ")"), round(summary(fit.nocc)$coefficients[, c(2, 5)], 2), summary(fit.nocc)$n)

        write.table(rbind(Result, Result.cc, Result.nocc), file = "NEJGenoGroup.xls", sep = "\t", quote = F, col.names = NA)



Result <- matrix(0, ncol = 9, nrow = 4, dimnames = list(c("AllUni", "AllAdj", "CCUni", "CCAdj"), c("lower.95", "upper.95", "coef", "exp_coef", "se_coef", "z", "P", "n", "Event")))

fit <- coxph(Surv(GL_time, GL_event) ~ NEJRiskHeta, data = GOCAR)

Result["AllUni", ] <- c(summary(fit)$conf.int["NEJRiskHeta", 3:4], summary(fit)$coefficients["NEJRiskHeta", ], summary(fit)$n, summary(fit)$nevent)

fit <- coxph(Surv(GL_time, GL_event) ~ NEJRiskHeta + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum, data = GOCAR)

Result["AllAdj", ] <- c(summary(fit)$conf.int["NEJRiskHeta", 3:4], summary(fit)$coefficients["NEJRiskHeta", ], summary(fit)$n, summary(fit)$nevent)

GOCARCC <- GOCAR[GOCAR$R2R == "Caucasian > Caucasian", ]

fit <- coxph(Surv(GL_time, GL_event) ~ NEJRiskHeta, data = GOCARCC)

Result["CCUni", ] <- c(summary(fit)$conf.int["NEJRiskHeta", 3:4], summary(fit)$coefficients["NEJRiskHeta", ], summary(fit)$n, summary(fit)$nevent)

fit <- coxph(Surv(GL_time, GL_event) ~ NEJRiskHeta + noLIMS1 + Induction2 + Donor_Status2 + HLA.mm.sum, data = GOCARCC)

Result["CCAdj", ] <- c(summary(fit)$conf.int["NEJRiskHeta", 3:4], summary(fit)$coefficients["NEJRiskHeta", ], summary(fit)$n, summary(fit)$nevent)

Tab <- data.frame(rbind(NA, Result))

rownames(Tab)[1] <- "Title"

Tab$Label <- factor(rownames(Tab), levels = rev(rownames(Tab)))

Tab$colour <- rep(c("white", "gray95"), 3)[1:nrow(Tab)]


Tab.man <- Tab

Tab.point <- Tab[which(Tab$upper.95 > 5), ]

Tab.point$upper.95 <- 5

Tab.man$upper.95[which(Tab.man$upper.95 > 5)] <- 5

fp <- ggplot(data=Tab.man, aes(x=Label, y=exp_coef, ymin=lower.95, ymax=upper.95)) +

        geom_vline(aes(xintercept = Label, colour = colour), size = 10) +

        geom_pointrange(color = Colors["Yellow", ]) +

        geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip

        geom_point(data = Tab.point, aes(x = Label, y = upper.95), shape = 62, color = Colors["Yellow", ], size = 4) +

	ylim(NA, 5) +

        scale_x_discrete(breaks = rownames(Tab.man), labels = c("", "All univariate", "All adjusted", "E->E univariate", "E->E adjusted")) +

        coord_flip() +  # flip coordinates (puts labels on y axis)

        xlab("rs893403 mismatch") + ylab("Hazard ratio (95% CI)") +

        theme_classic() + scale_colour_identity() + # use a white background

        theme(axis.ticks = element_blank(), axis.text.y = element_text(size = 10))


Tab$Range <- paste0(round(Tab$exp_coef, 2), " (", round(Tab$lower.95, 2), ", ", round(Tab$upper.95, 2), ")")

Tab$Range[Tab$Range == "NA (NA, NA)"] <- NA

Tab$P <- as.character(round(Tab$P, 3))

Tab$P[Tab$P == "0"] <- "<0.001"

Tab$P[1] <- "P value"

Tab$Range[1] <- "Hazard ratio (95% CI)"

fp2 <- ggplot(data = Tab, aes(y = Label)) +

        geom_hline(aes(yintercept = Label, colour = colour), size = 10) +

        #geom_text(aes(x = 0, label = Label), hjust = 0) +

        geom_text(aes(x = 0, label = P), hjust = 0, size = 3) +

        geom_text(aes(x = 3, label = Range), hjust = 1, size = 3) +

        scale_colour_identity() +

        theme_void() +

        theme(plot.margin = margin(5, 0, 32, 0))


Files <- c("All_LIMS1noRS_multi_iterate.xls", "CC_LIMS1noRS_multi_iterate.xls", "Stratified_LIMS1_multi_iterate.All.xls", "Stratified_LIMS1_multi_iterate.CC.xls")

Tabs <- lapply(Files, function(x) read.delim(file = x, row.names = 1))

Tabs <- lapply(Tabs, function(x) {x$FDR.BH = p.adjust(x$Pr...z...snp, method = "BH"); x})

names(Tabs) <- c("AdjustAll", "AdjustC2C", "StratAll", "StratC2C")

Tabs <- lapply(Tabs, function(x) x[(x$Pr...z...snp <= 0.01), ])

GeneS <- lapply(Tabs, function(x) rownames(x))


pp3 <- ggVennDiagram(GeneS, label_alpha = 0, category.names = c("All \nadjusted","E->E \nadjusted","All \nstratified", "E->E \nstratified"), label = "count", label_size = 4, set_size = 4, color = "black") + theme(legend.position = "null", plot.margin = margin(10, 10, 10, 20)) + scale_fill_gradient(low=brewer.pal(9, "YlOrRd")[2],high = brewer.pal(9, "YlOrRd")[6])


Overlap <- Reduce(intersect, GeneS)

for(nn in names(Tabs)){

        names(Tabs[[nn]]) <- paste0(nn, "_", names(Tabs[[nn]]))

}

Tmp <- Reduce(cbind, lapply(Tabs, function(x) x[Overlap, ]))

LIMS1SNPs <- read.delim(file = "LIMS1.SNPs", row.names = 1, as.is = T)

Tmp <- cbind(Tmp, LIMS1SNPs[rownames(Tmp), c("Chr", "Start", "Ref", "Alt", "Func.refGene")])

write.table(Tmp, file = "IterateOverlap_P01.xls", sep = "\t", quote = F, col.names = NA)

#ggsave("Test.pdf", pp3, width = 9)

Tab <- read.delim(file = "All_LIMS1noRS_multi_iterate.xls", row.names = 1, as.is = T)

Ann <- read.delim(file = "Variant_Matrix.xls.LIMS1.comm.heta.ann", row.names = 1, as.is = T, header = F)

Tab <- Tab[rownames(Ann), ]

Tab$Chr <- Ann[rownames(Tab), 1]

Tab$Coord <- Ann[rownames(Tab), 2]

Tab$Ref <- Ann[rownames(Tab), 4]

Tab$Alt <- Ann[rownames(Tab), 5]

Tab$Group <- Ann[rownames(Tab), 6]

Tab <- Tab[order(Tab$Coord, decreasing = F), ]

Tab$Name <- factor(rownames(Tab), rev(rownames(Tab)))

pp4 <- ggplot(Tab, aes(x = Name, y = -log10(Pr...z...snp), fill = Group)) +

	geom_bar(stat = "identity") +

	scale_fill_manual(values = Colors[c("Red", "Blue", "Yellow"), ]) +

	scale_y_continuous(expand = c(0, 0)) +

	xlab("") + ylab("-log10(P)") +

	theme_classic() +

	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank())

ForLocusZoom <- data.frame(Chromosome = Tab$Chr, Position = Tab$Coord, Ref.allele = Tab$Ref, Alt.allele = Tab$Alt, p.value = Tab$Pr...z...snp)
	
write.table(ForLocusZoom, file = "ForLocusZoom.tab", sep = "\t", row.names = F, quote = F)

#ggsave("FigNEJ.pdf", ggarrange(ggarrange(ggarrange(pp2$plot, pp2$table, nrow = 2, heights = c(2, 1)), pp3, widths = c(2, 2.5), labels = c("A", "C")), ggarrange(fp, fp2, ncol = 2, widths = c(3, 1), labels = "B"), ggarrange(pp4, nrow = 1, ncol = 1, labels = "D"), nrow = 3, heights = c(2, 1, 1)), width = 11, height = 10)

ggsave("FigNEJ.pdf", ggarrange(ggarrange(ggarrange(pp2$plot, pp2$table, nrow = 2, heights = c(2, 1)), pp3, widths = c(2, 2.5), labels = c("A", "C")), ggarrange(fp, fp2, ncol = 2, widths = c(3, 1), labels = "B"), ggarrange(ggplot() + theme_void(), nrow = 1, ncol = 1, labels = "D"), nrow = 3, heights = c(2, 1, 1)), width = 11, height = 10)


