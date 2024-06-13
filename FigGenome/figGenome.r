
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(survival)

scaleFUN2 <- function(x) sprintf("%.2f", x)

scaleFUN1 <- function(x) sprintf("%.1f", x)


GOCAR <- read.delim(file = "GOCAR.xls", row.names = 1, as.is = T)

CTOT <- read.delim(file = "CTOT.xls", row.names = 1, as.is = T)

Colors <- read.delim(file = "ColorLibrary", row.names = 1, as.is = T, header = F)

p <- ggplot(GOCAR, aes(x=All_freq)) + 

  	geom_density(fill = brewer.pal(8, "Set2")[1]) +

	theme_classic()

p2 <- ggplot(CTOT, aes(x=All_freq)) +

        geom_density(fill = brewer.pal(8, "Set2")[1]) +

        theme_classic()


ggsave("GOCAR_All_freq.density.pdf", arrangeGrob(grobs = list(p, p2), ncol = 2), height = 4, width = 8)


Race <- read.table(file = "Race_Genetic.xls", header = FALSE, sep = "\t", as.is = T)

names(Race) <- c("R2R", "Count", "Simplify", "Category")

Race$Simplify <- gsub("B", "AA", Race$Simplify)

GOCAR$R2R <- paste(GOCAR$Genetic_Donor_Race, GOCAR$Genetic_Rec_Race, sep=" > ")

GOCAR$Category <- Race$Category[match(GOCAR$R2R, Race$R2R)]

GOCAR$Simplify <- Race$Simplify[match(GOCAR$R2R, Race$R2R)]


GOCAR$Category <- factor(GOCAR$Category, levels = c("Inter", "Intra"))

levels(GOCAR$Category) <- c("Inter ancestry", "Intra ancestry")

rownames(Race) <- Race$Simplify

Race$Median <- tapply(GOCAR$All, GOCAR$Simplify, median)[rownames(Race)]

Race <- Race[with(Race, order(Category, Median)), ]

GOCAR$Simplify <- factor(GOCAR$Simplify, levels = Race$Simplify)


p1 <- ggplot(GOCAR, aes(x=Simplify, y=All, fill = Category)) +

        geom_boxplot(outlier.shape = NA) +

        geom_jitter(shape=20, position=position_jitter(0.2)) +

        theme_classic() + scale_y_continuous(expand = c(0, 0), labels = scaleFUN1) +

        scale_fill_manual(values = Colors[c("Red", "Blue"), ]) +

        theme(axis.title = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust=1, face = "bold"), legend.title = element_blank(), plot.margin = margin(10, 20, 0, 10, "pt")) +

	#theme(axis.title = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust=1, face = "bold"), legend.title = element_blank()) +

        xlab("") +

        ylab("Normalized mismatch score") +

        labs(fill = "Race")


Cor <- cor.test(GOCAR$PI_HAT, GOCAR$All)

p2 <- ggplot(GOCAR, aes(x=PI_HAT, y=All, color = Category)) +

        geom_point(size = 2, shape = 20) +

        geom_smooth(method=lm, linetype="dashed", color = Colors["SteelDark", ]) +

        scale_color_manual(values = Colors[c("Red", "Blue"), ]) +

        annotate("text", x = Inf, y = Inf, label = paste0("P < 0.001", "\n", "R = ", round(Cor$estimate, 3)), vjust = 1, hjust = 1, size = 3) +

        theme_classic() + scale_y_continuous(expand = c(0, 0), labels = scaleFUN1) + scale_x_continuous(expand = c(0, 0), labels = scaleFUN2, limits = c(0, 1.02)) +

        ylab("Normalized mismatch score") +

        xlab("pIBD") +

        theme(legend.title = element_blank(), axis.title = element_text(size = 10), legend.position = "none", plot.margin = margin(10, 10, 0, 20, "pt"))


write.table(GOCAR[, c("PI_HAT", "All", "Simplify", "Category", "Exon2", "noExon2", "AllWithoutNoSy", "NoSyAllWithoutTrmem", "trmem")], file = "FigGenomeGocar.xls", sep = "\t", quote = F, col.names = NA)


CTOT$Category <- sapply(CTOT$R2R, function(x) {tmp <- unlist(strsplit(as.character(x), split = " > ")); if(tmp[1] == tmp[2]){"Intra"} else {"Inter"}})

CTOT$Category <- factor(CTOT$Category, levels = c("Inter", "Intra"))

levels(CTOT$Category) <- c("Inter ancestry", "Intra ancestry")

Race <- read.delim(file = "Race.xls", as.is = T)

rownames(Race) <- Race$Simplify

CTOT$Simplify <- Race$Simplify[match(CTOT$R2R, Race$R2R)]

Race$Median <- tapply(CTOT$All, CTOT$Simplify, median)[rownames(Race)]

Race <- Race[with(Race, order(Category, Median)), ]

CTOT$Simplify <- factor(CTOT$Simplify, levels = Race$Simplify)


p3 <- ggplot(CTOT, aes(x=Simplify, y=All, fill = Category)) +

        geom_boxplot(outlier.shape = NA) +

        geom_jitter(shape=20, position=position_jitter(0.2)) +

        theme_classic() + scale_y_continuous(labels = scaleFUN1) +

        scale_fill_manual(values = Colors[c("Red", "Blue"), ]) +

        theme(axis.title = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust=1, face = "bold"), legend.title = element_blank(), plot.margin = margin(10, 20, 0, 10, "pt")) +

        xlab("") +

        ylab("Normalized mismatch score") +

        labs(fill = "Race")


Cor <- cor.test(CTOT$pIBD, CTOT$All)

p4 <- ggplot(CTOT, aes(x=pIBD, y=All, color = Category)) +

        geom_point(size = 2, shape = 20) +

        geom_smooth(method=lm, linetype="dashed", color = Colors["SteelDark", ]) +

        #scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +

        scale_color_manual(values = Colors[c("Red", "Blue"), ]) +

        annotate("text", x = Inf, y = Inf, label = paste0("P < 0.001", "\n", "R = ", round(Cor$estimate, 3)), vjust = 1, hjust = 1, size = 3) +

        theme_classic() + scale_x_continuous(expand = c(0, 0), labels = scaleFUN2) + scale_y_continuous(labels = scaleFUN1) +

        ylab("Normalized mismatch score") +

        xlab("pIBD") +

        theme(legend.title = element_blank(), axis.title = element_text(size = 10), legend.position = "none", plot.margin = margin(10, 10, 0, 20, "pt"))

write.table(CTOT[, c("pIBD", "All", "Simplify", "Category", "Exon2", "noExon2", "AllWithoutNoSy", "NoSyAllWithoutTrmem", "trmem")], file = "FigGenomeCTOT.xls", sep = "\t", quote = F, col.names = NA)


#ggsave("FigGenome.pdf", ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"), width = 8, height = 6)


        #Variates <- c("All", "Exon2", "noExon2", "AllWithoutNoSy", "NoSyAllWithoutTrmem", "trmem")

	#Variates <- c("trmem", "NoSyAllWithoutTrmem", "AllWithoutNoSy")

	Variates <- c("All", "Exon2", "noExon2")

        Result <- matrix(0, ncol = 9, nrow = length(Variates)*2, dimnames = list(c(paste0(Variates, "_uni"), paste0(Variates, "_adj")), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "Mean", "SD")))

        for(vv in Variates){

                fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~ ", vv)), data = GOCAR)

                Result[paste0(vv, "_uni"), ] <- c(summary(fit)$conf.int[vv, 3:4], summary(fit)$coefficients[vv, ], mean(GOCAR[[paste0(vv, "_freq")]]), sd(GOCAR[[paste0(vv, "_freq")]]))

                fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~ ", vv, " + Induction2 + Donor_Status2 + HLA.mm.sum")), data = GOCAR)

                Result[paste0(vv, "_adj"), ] <- c(summary(fit)$conf.int[vv, 3:4], summary(fit)$coefficients[vv, ], mean(GOCAR[[paste0(vv, "_freq")]]), sd(GOCAR[[paste0(vv, "_freq")]]))

        }


        Result2 <- matrix(0, ncol = 9, nrow = length(Variates)*2, dimnames = list(c(paste0(Variates, "_uni"), paste0(Variates, "_adj")), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "Mean", "SD")))

        for(vv in Variates){

                fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~ ", vv)), data = CTOT)

                Result2[paste0(vv, "_uni"), ] <- c(summary(fit)$conf.int[vv, 3:4], summary(fit)$coefficients[vv, ], mean(CTOT[[paste0(vv, "_freq")]]), sd(CTOT[[paste0(vv, "_freq")]]))

                fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~ ", vv, " + induct + DTYPE + HLAmm.6variable")), data = CTOT)

                Result2[paste0(vv, "_adj"), ] <- c(summary(fit)$conf.int[vv, 3:4], summary(fit)$coefficients[vv, ], mean(CTOT[[paste0(vv, "_freq")]]), sd(CTOT[[paste0(vv, "_freq")]]))

        }

#stop()

#rownames(Result) <- paste(rownames(Result), "GOCAR", sep = "_")

#rownames(Result2) <- paste(rownames(Result2), "CTOT01", sep = "_")

Tab <- data.frame(rbind(NA, Result[1:3, ], NA, Result[4:6, ], NA, Result2[1:3, ], NA, Result2[4:6, ]))

names(Tab) <- c("lower.95", "upper.95", "coef", "exp_coef", "se_coef", "z", "P", "Mean", "SD")

Tab$Label <- c("GOCAR univariate", "All SNPs univariate GOCAR", "Exonic SNPs univariate GOCAR", "Non-exonic SNPs univariate GOCAR", "GOCAR adjusted", "All SNPs adjusted GOCAR", "Exonic SNPs adjusted GOCAR", "Non-exonic SNPs adjusted GOCAR", "CTOT univariate", "All SNPs univariate CTOT", "Exonic SNPs univariate CTOT", "Non-exonic SNPs univariate CTOT", "CTOT adjusted", "All SNPs adjusted CTOT", "Exonic SNPs adjusted CTOT", "Non-exonic SNPs adjusted CTOT")

Tab$Label <- factor(Tab$Label, levels = rev(Tab$Label))

Tab$colour <- rep(c("white", "gray95"), 8)

Tab.man <- Tab

Tab.point <- Tab[which(Tab$upper.95 > 3.5), ]

Tab.point$upper.95 <- 3.5

Tab.man$upper.95[which(Tab.man$upper.95 > 3.5)] <- 3.5

#stop()

fp <- ggplot(data=Tab.man, aes(x=Label, y=exp_coef, ymin=lower.95, ymax=upper.95)) +

	geom_vline(aes(xintercept = Label, colour = colour), size = 6) +

        geom_pointrange(color = Colors["Yellow", ]) + 

        geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip

	geom_point(data = Tab.point, aes(x = Label, y = upper.95), shape = 62, color = Colors["Yellow", ], size = 4) +

	ylim(NA, 3.5) +

	scale_x_discrete(breaks=c("GOCAR univariate", "All SNPs univariate GOCAR", "Exonic SNPs univariate GOCAR", "Non-exonic SNPs univariate GOCAR", "GOCAR adjusted", "All SNPs adjusted GOCAR", "Exonic SNPs adjusted GOCAR", "Non-exonic SNPs adjusted GOCAR", "CTOT univariate", "All SNPs univariate CTOT", "Exonic SNPs univariate CTOT", "Non-exonic SNPs univariate CTOT", "CTOT adjusted", "All SNPs adjusted CTOT", "Exonic SNPs adjusted CTOT", "Non-exonic SNPs adjusted CTOT"), labels=c("GOCAR univariate", "All SNPs", "Exonic SNPs", "Non-exonic SNPs", "GOCAR adjusted", "All SNPs", "Exonic SNPs", "Non-exonic SNPs", "CTOT univariate", "All SNPs", "Exonic SNPs", "Non-exonic SNPs", "CTOT adjusted", "All SNPs", "Exonic SNPs", "Non-exonic SNPs")) +

        coord_flip() +  # flip coordinates (puts labels on y axis)

        xlab("") + ylab("Hazard ratio (95% CI)") +

        theme_classic() + scale_colour_identity() + # use a white background

	theme(axis.ticks = element_blank())


Tab$Range <- paste0(round(Tab$exp_coef, 2), " (", round(Tab$lower.95, 2), ", ", round(Tab$upper.95, 2), ")")

Tab$Range[Tab$Range == "NA (NA, NA)"] <- NA

Tab$P <- as.character(round(Tab$P, 3))

Tab$P[Tab$P == "0"] <- "<0.001"

Tab[1, c("P", "Range")] <- c("P value", "Hazard ratio (95% CI)")

fp2 <- ggplot(data = Tab, aes(y = Label)) +

        geom_hline(aes(yintercept = Label, colour = colour), size = 6) +

  	#geom_text(aes(x = 0, label = Label), hjust = 0) +

  	geom_text(aes(x = 0, label = P), hjust = 0, size = 3) +

 	geom_text(aes(x = 3, label = Range), hjust = 1, size = 3) +

  	scale_colour_identity() +

  	theme_void() +

  	theme(plot.margin = margin(5, 0, 32, 0))

#ggsave("Test.pdf", ggarrange(fp, fp2, ncol = 2, widths = c(3, 1)), width = 12)

#stop()

#ggsave("FigGenome.pdf", ggarrange(ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"), ggarrange(fp, ggplot() + theme_void(), nrow = 1, ncol = 2, labels = "E", widths = c(2, 1)), nrow = 2, heights = c(2, 1)), width = 8, height = 10)

ggsave("FigGenome.pdf", ggarrange(ggarrange(p1, ggarrange(p2, nrow = 2, heights = c(7, 1)), p3, ggarrange(p4, nrow = 2, heights = c(7, 1)), labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"), ggarrange(ggarrange(fp, fp2, ncol = 2, widths = c(3, 1)), nrow = 1, labels = "E"), nrow = 2, heights = c(2, 1)), width = 9, height = 10)

write.table(Tab, file = "FigGenomeForrest.xls", sep = "\t", quote = F, col.names = NA)

#stop()

	GOCARCC <- GOCAR[GOCAR$R2R == "Caucasian > Caucasian", ]

	CTOTCC <- CTOT[CTOT$R2R == "EUR > EUR", ]

        Variates <- c("All", "Exon2", "noExon2")

        Result <- matrix(0, ncol = 9, nrow = length(Variates)*2, dimnames = list(c(paste0(Variates, "_uni"), paste0(Variates, "_adj")), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "Mean", "SD")))

        for(vv in Variates){

                fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~ ", vv)), data = GOCARCC)

                Result[paste0(vv, "_uni"), ] <- c(summary(fit)$conf.int[vv, 3:4], summary(fit)$coefficients[vv, ], mean(GOCARCC[[paste0(vv, "_freq")]]), sd(GOCARCC[[paste0(vv, "_freq")]]))

                fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~ ", vv, " + Induction2 + Donor_Status2 + HLA.mm.sum")), data = GOCARCC)

                Result[paste0(vv, "_adj"), ] <- c(summary(fit)$conf.int[vv, 3:4], summary(fit)$coefficients[vv, ], mean(GOCARCC[[paste0(vv, "_freq")]]), sd(GOCARCC[[paste0(vv, "_freq")]]))

        }

        Result2 <- matrix(0, ncol = 9, nrow = length(Variates)*2, dimnames = list(c(paste0(Variates, "_uni"), paste0(Variates, "_adj")), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "Mean", "SD")))

        for(vv in Variates){

                fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~ ", vv)), data = CTOTCC)

                Result2[paste0(vv, "_uni"), ] <- c(summary(fit)$conf.int[vv, 3:4], summary(fit)$coefficients[vv, ], mean(CTOTCC[[paste0(vv, "_freq")]]), sd(CTOTCC[[paste0(vv, "_freq")]]))

                fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~ ", vv, " + induct + DTYPE + HLAmm.6variable")), data = CTOTCC)

                Result2[paste0(vv, "_adj"), ] <- c(summary(fit)$conf.int[vv, 3:4], summary(fit)$coefficients[vv, ], mean(CTOTCC[[paste0(vv, "_freq")]]), sd(CTOTCC[[paste0(vv, "_freq")]]))

        }

rownames(Result) <- paste(rownames(Result), "GOCAR", sep = "_")

rownames(Result2) <- paste(rownames(Result2), "CTOT01", sep = "_")

#Tab <- data.frame(rbind(Result, Result2))

Tab <- data.frame(rbind(NA, Result[1:3, ], NA, Result[4:6, ], NA, Result2[1:3, ], NA, Result2[4:6, ]))

names(Tab) <- c("lower.95", "upper.95", "coef", "exp_coef", "se_coef", "z", "P", "Mean", "SD")

Tab$Label <- c("GOCAR univariate", "All SNPs univariate GOCAR", "Exonic SNPs univariate GOCAR", "Non-exonic SNPs univariate GOCAR", "GOCAR adjusted", "All SNPs adjusted GOCAR", "Exonic SNPs adjusted GOCAR", "Non-exonic SNPs adjusted GOCAR", "CTOT univariate", "All SNPs univariate CTOT", "Exonic SNPs univariate CTOT", "Non-exonic SNPs univariate CTOT", "CTOT adjusted", "All SNPs adjusted CTOT", "Exonic SNPs adjusted CTOT", "Non-exonic SNPs adjusted CTOT")

Tab$Label <- factor(Tab$Label, levels = rev(Tab$Label))

Tab$colour <- rep(c("white", "gray95"), 8)

Tab.man <- Tab

Tab.point <- Tab[which(Tab$upper.95 > 3.5), ]

Tab.point$upper.95 <- 3.5

Tab.man$upper.95[which(Tab.man$upper.95 > 3.5)] <- 3.5

fp <- ggplot(data=Tab.man, aes(x=Label, y=exp_coef, ymin=lower.95, ymax=upper.95)) +

        geom_vline(aes(xintercept = Label, colour = colour), size = 6) +

        geom_pointrange(color = Colors["Yellow", ]) +

        geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip

        geom_point(data = Tab.point, aes(x = Label, y = upper.95), shape = 62, color = Colors["Yellow", ], size = 4) +

        ylim(NA, 3.5) +

        scale_x_discrete(breaks=c("GOCAR univariate", "All SNPs univariate GOCAR", "Exonic SNPs univariate GOCAR", "Non-exonic SNPs univariate GOCAR", "GOCAR adjusted", "All SNPs adjusted GOCAR", "Exonic SNPs adjusted GOCAR", "Non-exonic SNPs adjusted GOCAR", "CTOT univariate", "All SNPs univariate CTOT", "Exonic SNPs univariate CTOT", "Non-exonic SNPs univariate CTOT", "CTOT adjusted", "All SNPs adjusted CTOT", "Exonic SNPs adjusted CTOT", "Non-exonic SNPs adjusted CTOT"), labels=c("GOCAR univariate", "All SNPs", "Exonic SNPs", "Non-exonic SNPs", "GOCAR adjusted", "All SNPs", "Exonic SNPs", "Non-exonic SNPs", "CTOT univariate", "All SNPs", "Exonic SNPs", "Non-exonic SNPs", "CTOT adjusted", "All SNPs", "Exonic SNPs", "Non-exonic SNPs")) +

        coord_flip() +  # flip coordinates (puts labels on y axis)

        xlab("") + ylab("Hazard ratio (95% CI)") +

        theme_classic() + scale_colour_identity() + # use a white background

        theme(axis.ticks = element_blank())

Tab$Range <- paste0(round(Tab$exp_coef, 2), " (", round(Tab$lower.95, 2), ", ", round(Tab$upper.95, 2), ")")

Tab$Range[Tab$Range == "NA (NA, NA)"] <- NA

Tab$P <- as.character(round(Tab$P, 3))

Tab$P[Tab$P == "0"] <- "<0.001"

Tab[1, c("P", "Range")] <- c("P value", "Hazard ratio (95% CI)")

fp2 <- ggplot(data = Tab, aes(y = Label)) +

        geom_hline(aes(yintercept = Label, colour = colour), size = 6) +

        #geom_text(aes(x = 0, label = Label), hjust = 0) +

        geom_text(aes(x = 0, label = P), hjust = 0, size = 3) +

        geom_text(aes(x = 3, label = Range), hjust = 1, size = 3) +

        scale_colour_identity() +

        theme_void() +

        theme(plot.margin = margin(5, 0, 32, 0))

ggsave("FigGenome.CC.pdf", ggarrange(fp, fp2, ncol = 2, widths = c(3, 1)), height = 4, width = 8)

write.table(Tab, file = "FigGenomeForrest.EE.xls", sep = "\t", quote = F, col.names = NA)


#stop()

DatList <- list(GOCAR = GOCAR, GOCARCC = GOCARCC, CTOT = CTOT, CTOTCC = CTOT)

scaleFUN <- function(x) sprintf("%.1f", x)

Variates <- c("All", "Exon2", "noExon2", "AllWithoutNoSy", "NoSyAllWithoutTrmem", "trmem")

Plots <- list()

for(dd in names(DatList)){

for(nn in Variates){

	for(nn2 in Variates){

		Cor <- cor.test(DatList[[dd]][[nn]], DatList[[dd]][[nn2]])

		pp <- ggplot(data = DatList[[dd]], aes_string(x = nn2, y = nn)) +

			geom_point(size = 2, shape = 20, color = Colors["Blue", ]) +

			#xlim(0, NA) + ylim(0, NA) +

        		geom_smooth(method=lm, linetype="dashed", color = Colors["SteelDark", ]) +

			annotate("text", x = 0, y = Inf, label = paste0("P = ", round(Cor$p.value, 3), "\n", "R = ", round(Cor$estimate, 3)), vjust = 1, hjust = 0, size = 3) +

        		theme_classic() + xlab("") + ylab("") +

			scale_y_continuous(labels=scaleFUN, limits = c(0, NA)) + scale_x_continuous(labels=scaleFUN, limits = c(0, NA)) +

        		theme(legend.title = element_blank(), axis.title = element_text(size = 10))

		Plots[[paste(nn, nn2, sep = "_")]] <- pp

	}

}

for(i in 2:length(Variates)){

	for(j in 1:(i - 1)){

		Plots[[length(Variates) * i - j + 1]] <- ggplot() + theme_void()

	}

}

Names <- c("All", "Exonic", "Non-exonic", "All-nsSNPs", "nsSNPs-trSNPs", "trSNPs")

Plots[1:length(Variates)] <- lapply(1:length(Variates), function(i) arrangeGrob(Plots[[i]], top=Names[i]))

Plots[((0:(length(Variates) - 1))*length(Variates) + 1)] <- lapply(1:length(Variates), function(i) arrangeGrob(Plots[[(i-1)*length(Variates) + 1]], left=Names[i]))

ggsave(paste0(dd, ".Cor_All.pdf"), arrangeGrob(grobs = Plots, ncol = 6), width = 12, height = 12)

}

