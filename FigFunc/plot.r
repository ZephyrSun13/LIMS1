
Args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(RColorBrewer)

Files <- unlist(read.table(file = "GO.lst", stringsAsFactors = FALSE))

Tabs <- lapply(Files, function(x) read.delim(file = x, stringsAsFactors = FALSE))

names(Tabs) <- sapply(Files, function(x) unlist(strsplit(x, split = "/"))[3])

Tabs <- lapply(Tabs, function(x) x[x$pvalue <= 0.1, ])

TabsGO <- Tabs


Files <- unlist(read.table(file = "KEGG.lst", stringsAsFactors = FALSE))

Tabs <- lapply(Files, function(x) read.delim(file = x, stringsAsFactors = FALSE))

names(Tabs) <- sapply(Files, function(x) unlist(strsplit(x, split = "/"))[3])

Tabs <- lapply(Tabs, function(x) x[x$pvalue <= 0.1, ])

TabsKEGG <- Tabs


Tabs <- list()

for(nn in names(TabsGO)){

	Tabs[[nn]] <- rbind(TabsGO[[nn]], TabsKEGG[[nn]])

}


GOPlot <- unlist(read.table(file = paste0(Args[1], "_GO.plot"), stringsAsFactors = FALSE, sep = "\t"))

for(nn in names(Tabs)){

	if(!any(Tabs[[nn]]$Description %in% GOPlot)){

		Tabs[[nn]] <- NULL

	}

}

Tab <- Reduce(rbind, lapply(names(Tabs), function(x) cbind(Tabs[[x]][Tabs[[x]]$Description %in% GOPlot, c("Description", "pvalue", "Count")], Group = x)))

#NameMap <- read.delim(file = "NameMap.xls", as.is = T, row.names = 1)

#NameMap <- NameMap[rownames(NameMap) %in% names(Tabs), , drop = F]

#Tab$Group <- factor(NameMap[as.character(Tab$Group), "Map"], levels = NameMap$Map)

#Tab$Description <- factor(as.character(Tab$Description), levels = rev(GOPlot))

Tab$Description <- factor(as.character(Tab$Description), levels = rev(GOPlot))

#Tab$Group <- factor(as.character(Tab$Group), c("Resident_Mac", "Mrc1_Mac", "Inflammatory_Mac", "IFN_Mac", "Proliferate_Mac", "Trem2_hi_Mac", "Infiltrate_Mac", "Ccr7_DC", "cDC1", "cDC2", "pDC"))

#Colors2 <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/ColorLibrary", as.is = T, row.names = 1, header = F)

pp <- ggplot(Tab, aes(x = Group, y = Description, color = -log(pvalue))) +

	geom_point(aes(size = Count)) +

	scale_size(range = c(5, 10)) +

	#scale_color_gradient(high = Colors2["SteelDark", ], low = Colors2["SteelLight", ]) +

        scale_color_gradient(high = brewer.pal(9, "YlOrRd")[9], low = brewer.pal(9, "YlOrRd")[4]) +

	theme_minimal() +

	xlab("") +

	ylab("") +

	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"), axis.text.y = element_text(size = 20, face = "bold"), legend.text = element_text(size = 15), legend.title = element_text(size = 15))

ggsave(file = paste0(Args[1], "_GOPlot.pdf"), pp, height = 12, width = 8)

