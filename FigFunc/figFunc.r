
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(ggpubr)

Colors <- read.delim(file = "ColorLibrary", row.names = 1, as.is = T, header = F)

Files <- unlist(read.table(file = "Exprs", as.is = TRUE))

Group <- read.delim(file = "Haplo.xls", row.names = 1, as.is = TRUE)

Group[Group == "0|1"] <- "1|0"

Group$rs893403 <- factor(Group$rs893403)

levels(Group$rs893403) <- c("G/G", "G/A", "A/A")

Group$rs893403 <- as.character(Group$rs893403)

Group$SNP30 <- gsub("\\|", "/", Group$SNP30)


rownames(Group) <- paste0("DICE_", rownames(Group))


names(Files) <- gsub(".csv.gz", "", sapply(Files, function(x) unlist(strsplit(x, split = "/"))[10]))

#write.table(names(Files), file = "NameMap.xls", sep = "\t", quote = F, row.names = F, col.names = F)

Plots1 <- list()

Plots2 <- list()

Tabs <- list()

Cells <- c("TREG_NAIVE", "CD4_NAIVE", "CD8_NAIVE", "B_NAIVE", "CD4_N_STIM", "CD8_N_STIM", "CLASSICAL_MONOCYTES", "NK_CD16POS", "NONCLASSICAL_MONOCYTES", "TFH", "TH1-17", "TH17", "TH1", "TH2", "TREG_MEMORY")

Cells <- c(Cells, setdiff(names(Files), Cells))

Genes <- c("GCC2", "LIMS1")

for(nn in Cells){

	Expr <- read.delim(file = Files[nn], row.names = 1, as.is = T)

	rownames(Expr) <- sapply(rownames(Expr), function(x) unlist(strsplit(x, split = "\\."))[1])

	Tab <- data.frame(LIMS1 = unlist(Expr["ENSG00000169756", ]), GCC2 = unlist(Expr["ENSG00000135968", ]), GCC2_AS1 = unlist(Expr["ENSG00000214184", ]), Haplotype = Group[names(Expr), "SNP30"], rs893403 = Group[names(Expr), "rs893403"])

	Tab$rs893403 <- factor(as.character(Tab$rs893403), levels = c("G/G", "G/A", "A/A"))

	for(gg in Genes){

		for(cc in c("Haplotype", "rs893403")){

			pp <- ggplot(Tab, aes_string(y = gg, x = cc, fill = cc)) +

				geom_boxplot(outlier.shape = NA) + 

				stat_compare_means(method = "anova") +

				geom_jitter(shape=16, position=position_jitter(0.2)) +

				ggtitle(nn) +

				scale_fill_manual(values = Colors[c("Blue", "Yellow", "Red"), ]) +

				theme_classic() +

				theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

			if(cc == "Haplotype"){

				Plots1[[paste0(nn, "_", gg, "_", cc)]] <- pp

			}else if(cc == "rs893403"){

				Plots2[[paste0(nn, "_", gg, "_", cc)]] <- pp

			}

		}

	}

	Tabs[[nn]] <- Tab

}

#saveRDS(Plots, file = "Plots.rds")

#ggsave("FigFunc.pdf", ggarrange(Plots2[[1]], Plots2[[2]], Plots2[[3]], Plots1[[1]], Plots1[[2]], Plots1[[3]], ncol = 6, nrow = 1), width = 15, height = 3)


Tabs <- lapply(names(Tabs), function(x) cbind(Tabs[[x]], Cell = x))

Tab <- Reduce(rbind, Tabs)


Colors <- c(brewer.pal(8, "Set2"), brewer.pal(9, "Set1"), brewer.pal(12, "Set3"))

#Colors <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set2"), brewer.pal(9, "Set1"))

Plots <- list()

Genes <- c("LIMS1", "GCC2")

NameMap <- unlist(read.table(file = "NameMap.xls", header = F, as.is = T))

Tab$Cell <- factor(as.character(Tab$Cell), levels = NameMap)

for(gg in Genes){

	pp <- ggplot(Tab, aes_string(x = "Cell", y = gg, fill = "Cell")) +

		geom_violin() +

		#geom_boxplot(outlier.shape = NA, width = 0.1) +

		scale_fill_manual(values = Colors) +

		theme_classic() +

		theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), legend.position = "none", axis.title.x = element_blank())

	Plots[[gg]] <- pp

}

write.table(Tab, file = "DICEeQTL.table.xls", quote = F, sep = "\t", col.names = NA)

ggsave("FigFunc.pdf", ggarrange(ggarrange(Plots2[[1]], Plots2[[3]], Plots2[[5]], Plots2[[2]], Plots2[[4]], Plots2[[6]], Plots1[[1]], Plots1[[3]], Plots1[[5]], Plots1[[2]], Plots1[[4]], Plots1[[6]], ncol = 6, nrow = 2), ggarrange(ggarrange(Plots[["GCC2"]], nrow = 1, ncol = 1), ggplot() + theme_void(), ncol = 2, widths = c(2, 1)), nrow = 2, ncol = 1, labels = c("A", "B"), heights = c(2, 1.5)), width = 15, height = 9)

ggsave("FigFuncSupp1.pdf", arrangeGrob(grobs = Plots1[7:length(Plots1)], ncol = 4), width = 12, height = 18)

ggsave("FigFuncSupp2.pdf", arrangeGrob(grobs = Plots2[7:length(Plots1)], ncol = 4), width = 12, height = 18)

ggsave("FigFuncSupp3.pdf", Plots[["LIMS1"]], height = 3)
