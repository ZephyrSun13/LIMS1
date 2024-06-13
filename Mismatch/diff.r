
#!/usr/bin/R

#library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)

SNPs <- read.table(file = args[1], header = FALSE, sep = "\t", stringsAsFactors = FALSE)

#SNPs <- read.table(file = "Sun.tped", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

SNPs <- SNPs[!(duplicated(SNPs[[2]])), ]

rownames(SNPs) <- SNPs[[2]]

SNPs.ann <- SNPs[, 1:4]

SNPs <- SNPs[, -(1:4)]

fam <- read.table(file = "/sc/orga/projects/GOCAR/Sun/5.DRWork/2.impute/all.tfam.gz", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

#fam <- cbind(fam, sapply(fam[[2]], function(x) unlist(strsplit(x, "___"))[1]))

names(SNPs) <- fam[[2]]

Info <- read.table(file = "/sc/orga/projects/GOCAR/Sun/5.DRWork/1.SNP/DRinfo.xls", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

Sams <- Info[, c("Rec_GT_ID", "Donor_GT_ID")]

#Sams <- sapply(unlist(Info[, c("Rec_GT_ID", "Donor_GT_ID")]), function(x) unlist(strsplit(x, "___"))[1])

#SNPs <- SNPs[, Sams]

Matr <- matrix(0, nrow = nrow(SNPs), ncol = nrow(Sams), dimnames = list(rownames(SNPs), rownames(Sams)))

count <- 0

for(i in rownames(SNPs)){

	count <- count + 1

	if(count %% 10000 == 0){print(count)}

	for(j in rownames(Sams)){

                if(SNPs[i, Sams[j, "Rec_GT_ID"]] == "0 0" || SNPs[i, Sams[j, "Donor_GT_ID"]] == "0 0"){

                        Matr[i, j] <- NA

                        next

                }

		Recs <- unique(unlist(strsplit(SNPs[i, Sams[j, "Rec_GT_ID"]], " ")))

		Dons <- unique(unlist(strsplit(SNPs[i, Sams[j, "Donor_GT_ID"]], " ")))

		if(length(setdiff(Dons, Recs)) > 0){

			if(length(Dons) == 1 && length(Recs) == 1){

				Matr[i, j] <- 2

			}else{

				Matr[i, j] <- 1
			

			}

		}

	}

}

#write.table(Matr, file = paste0"Variant_Matrix.xls", col.names = NA, quote = FALSE, sep = "\t")

write.table(Matr, file = paste0(args[1], "Variant_Matrix.xls"), col.names = NA, quote = FALSE, sep = "\t")

#pdf("Heatmap.pdf")
#
#	pheatmap(Matr, border_color = NA, cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, legend = FALSE)
#
#dev.off()

