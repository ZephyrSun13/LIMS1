
library(pheatmap)

Tab <- read.delim(file = "SNP30.tab", header = FALSE, row.names = 1, as.is = TRUE)

fam <- read.table(file = "/sc/arion/projects/GOCAR/Sun/5.DRWork/2.impute/all.tfam.gz", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

Sams <- gsub("Db", "DB", fam[[1]])

Sams <- sapply(Sams, function(x) unlist(strsplit(x, split = "___"))[1])

names(Tab) <- Sams

Sams <- read.table(file = "/sc/arion/projects/GOCAR/Sun/5.DRWork/2.impute/DRinfo.xls", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)


Tab <- Tab[, c(Sams$Recipient, Sams$Donor)]


tt <- apply(Tab, 2, function(x) length(unique(x)))

table(tt)

Haplo <- apply(Tab, 2, function(x) {tmp <- table(x); names(tmp[tmp==max(tmp)])})

TabTmp <- Tab

TabTmp[TabTmp=="./."] <- NA

TabTmp[TabTmp=="0/0"] <- 0

TabTmp[TabTmp=="0/1"] <- 1

TabTmp[TabTmp=="1/1"] <- 2

TabTmp <- apply(TabTmp, 2, function(x) as.numeric(x))

rownames(TabTmp) <- rownames(Tab)

pdf("Haplo.heat.GoCAR.pdf", width = 20)

        pheatmap(TabTmp, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = F, show_rownames = T)

dev.off()


TabTmp2 <- t(apply(Tab, 1, function(x) sapply(x, function(y) as.integer(unlist(strsplit(y, split = "/"))))))

write.table(TabTmp2, file = "Haplo.heat2.GoCAR.xls", sep = "\t", quote = F, col.names = NA)

#TabTmp2[TabTmp2=="./."] <- NA

pdf("Haplo.heat2.GoCAR.pdf", width = 20)

        pheatmap(TabTmp2, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = F, show_rownames = T)

dev.off()


Clinical <- read.delim(file = "GOCAR.xls", row.names = 1, as.is = T)

table(Clinical$Genetic_Rec_Race)

table(Clinical$Genetic_Donor_Race)


Races <- c(Clinical$Genetic_Rec_Race, Clinical$Genetic_Donor_Race)

names(Races) <- c(rownames(Clinical), gsub("Rb", "DB", rownames(Clinical)))


MAF <- t(apply(TabTmp, 1, function(x) tapply(x, Races[colnames(TabTmp)], function(y) sum(y, na.rm = T)/(2*length(y)))))

MAFGOCAR <- MAF


Fam <- unlist(read.table(file = "/sc/arion/projects/GOCAR/Sun/5.DRWork/4.CTOT2/3.Diff/Sam", stringsAsFactors = FALSE))

Fam <- gsub("\\.2", "", Fam)

Fam <- sapply(Fam, function(x) paste(unlist(strsplit(x, split = "_"))[-1], collapse = "_"))

Rec <- grep("Recipient", Fam, value = TRUE)

Rec <- data.frame(ID = gsub("_Recipient", "", Rec), Rec = Rec, stringsAsFactors = FALSE)

Don <- grep("Donor", Fam, value = TRUE)

Don <- data.frame(ID = gsub("_Donor", "", Don), Don = Don, stringsAsFactors = FALSE)

Sams <- merge(Rec, Don, by = 1)

rownames(Sams) <- Sams$ID

CTOT <- read.delim(file = "../../4.CTOT2/4.SurvivalHaplo/Clinical_Haplo.xls", row.names = 1, as.is = T)

Sams <- Sams[rownames(Sams) %in% rownames(CTOT), ]

Tab <- read.delim(file = "SNP30.CTOT.vcf.GT.FORMAT.withID", header = TRUE, row.names = 1, as.is = TRUE, check.names = FALSE)

Tab <- Tab[, -(1:2)]

names(Tab) <- Fam

Tab <- Tab[, names(Tab) %in% unlist(Sams)]


tt <- apply(Tab, 2, function(x) length(unique(x)))

table(tt)


TabTmp <- Tab

TabTmp[TabTmp=="0|0"] <- 0

TabTmp[TabTmp=="0|1"] <- 1

TabTmp[TabTmp=="1|0"] <- 1

TabTmp[TabTmp=="1|1"] <- 2

TabTmp <- apply(TabTmp, 2, function(x) as.numeric(x))

rownames(TabTmp) <- rownames(Tab)

pdf("Haplo.heat.CTOT.pdf", width = 20)

        pheatmap(TabTmp, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = F, show_rownames = T)

dev.off()


TabTmp2 <- t(apply(Tab, 1, function(x) sapply(x, function(y) as.integer(unlist(strsplit(y, split = "\\|"))))))

write.table(TabTmp2, file = "Haplo.heat2.CTOT.xls", sep = "\t", quote = F, col.names = NA)

#TabTmp2[TabTmp2=="./."] <- NA

pdf("Haplo.heat2.CTOT.pdf", width = 20)

        pheatmap(TabTmp2, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = F, show_rownames = T)

dev.off()


table(CTOT$R_G_Race)

table(CTOT$D_G_Race)


Races <- c(CTOT$R_G_Race, CTOT$D_G_Race)

names(Races) <- c(paste0(rownames(CTOT), "_Recipient"), paste0(rownames(CTOT), "_Donor"))

table(Races)


MAF <- t(apply(TabTmp, 1, function(x) tapply(x, Races[colnames(TabTmp)], function(y) sum(y, na.rm = T)/(2*length(y)))))

MAFCTOT <- MAF

MAFAll <- cbind(MAFGOCAR, MAFCTOT)

write.table(MAFAll, file = "MAFAll.xls", sep = "\t", quote = F, col.names = NA)


## 1000G

Tab <- read.delim(file = "SNP30.1000G.vcf.GT.FORMAT.withID", header = TRUE, as.is = TRUE, check.names = FALSE)

Tab <- Tab[-which(Tab$ID == "chr2:108591293")[1], ]

rownames(Tab) <- Tab$ID

Sams <- unlist(read.table(file = "Sam2504", as.is = T))

Tab <- Tab[, Sams]

TabTmp <- Tab

TabTmp[TabTmp=="0|0"] <- 0

TabTmp[TabTmp=="0|1"] <- 1

TabTmp[TabTmp=="1|0"] <- 1

TabTmp[TabTmp=="1|1"] <- 2


TabTmp <- apply(TabTmp, 2, function(x) as.numeric(x))

rownames(TabTmp) <- rownames(Tab)

pdf("Haplo.heat.1000G.pdf", width = 20)

        pheatmap(TabTmp, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = F, show_rownames = T)

dev.off()


TabTmp2 <- t(apply(Tab, 1, function(x) sapply(x, function(y) as.integer(unlist(strsplit(y, split = "\\|"))))))

write.table(TabTmp2, file = "Haplo.heat2.1000G.xls", sep = "\t", quote = F, col.names = NA)

#TabTmp2[TabTmp2=="./."] <- NA

pdf("Haplo.heat2.1000G.pdf", width = 20)

        pheatmap(TabTmp2, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = F, show_rownames = T)

dev.off()

