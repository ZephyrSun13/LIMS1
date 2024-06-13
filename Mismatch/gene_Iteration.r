
library("survival")
library("survminer")

args <- commandArgs(trailingOnly = TRUE)

Tab <- read.table(file = args[1], sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

Clinical <- read.table(file = "/sc/arion/projects/GOCAR/zhangz05/data_to_Zeguo/phenotype_588samples_RD_4admix_update030217.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

Clinical$Genetic_Rec_Race <- NA
Clinical$Genetic_Rec_Race[(!is.na(Clinical$Rec_pAFR)) & (!is.na(Clinical$Rec_pEUR))] <- "Hispanic"
Clinical$Genetic_Rec_Race[Clinical$Rec_pAFR >= 0.6] <- "African American"
Clinical$Genetic_Rec_Race[Clinical$Rec_pEUR >= 0.9] <- "Caucasian"
Clinical$Genetic_Rec_Race[Clinical$Rec_pASN_AMR >= 0.9] <- "Asian"

Clinical$Genetic_Donor_Race <- NA
Clinical$Genetic_Donor_Race[(!is.na(Clinical$Donor_pAFR)) & (!is.na(Clinical$Donor_pEUR))] <- "Hispanic"
Clinical$Genetic_Donor_Race[Clinical$Donor_pAFR >= 0.6] <- "African American"
Clinical$Genetic_Donor_Race[Clinical$Donor_pEUR >= 0.9] <- "Caucasian"
Clinical$Genetic_Donor_Race[Clinical$Donor_pASN_AMR >= 0.9] <- "Asian"


DR_gene_all <- read.table(file = "/sc/arion/projects/GOCAR/Sun/5.DRWork/2.impute/Diff/Variant_Matrix.xls.all.freq.freq.xls", row.names = 1, header = FALSE, sep = "\t")

rownames(DR_gene_all) <- gsub("X", "", rownames(DR_gene_all))

Result_all <- matrix(NA, nrow = nrow(Tab), ncol = 18, dimnames = list(rownames(Tab), c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "Mean0", "Mean1", "SD0", "SD1", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "Mean0", "Mean1", "SD0", "SD1")))

#stop()

for(nn in rownames(Tab)){

	if((which(rownames(Tab) == nn) %% 1000) == 0){

		print(which(rownames(Tab) == nn))

	}

	DR_gene <- data.frame(t(Tab[nn, ]))

	DR_nongene <- DR_gene_all[rownames(DR_gene), , drop = FALSE] - DR_gene

		if(IQR(DR_gene[rownames(DR_gene),]) != 0){
	
			Clinical_tt <- data.frame(Clinical[rownames(DR_gene), ], gene = DR_gene[rownames(DR_gene),]/IQR(DR_gene[rownames(DR_gene),]), nongene = DR_nongene[rownames(DR_gene),]/IQR(DR_nongene[rownames(DR_gene),]), R2R = paste(Clinical[rownames(DR_gene), "Genetic_Donor_Race"], Clinical[rownames(DR_gene), "Genetic_Rec_Race"], sep=" > "))

		}else{

                        Clinical_tt <- data.frame(Clinical[rownames(DR_gene), ], gene = DR_gene[rownames(DR_gene),], nongene = DR_nongene[rownames(DR_gene),]/IQR(DR_nongene[rownames(DR_gene),]), R2R = paste(Clinical[rownames(DR_gene), "Genetic_Donor_Race"], Clinical[rownames(DR_gene), "Genetic_Rec_Race"], sep=" > "))

		}
		
		fit <- coxph(Surv(graft_loss_time_update, death_censored_graft_loss_event) ~ gene + nongene, data = Clinical_tt)

		Result <- summary(fit)$coefficients[c("gene", "nongene"), ]

		#Result <- cbind(Result, Mean = c(mean(DR_gene[[1]]), mean(DR_nongene[[1]])), SD = c(sd(DR_gene[[1]]), sd(DR_nongene[[1]])))

		Result <- cbind(Result, Mean = rbind(tapply(DR_gene[rownames(Clinical_tt), ], Clinical_tt$death_censored_graft_loss_event, mean), tapply(DR_nongene[rownames(Clinical_tt), ], Clinical_tt$death_censored_graft_loss_event, mean)), SD = rbind(tapply(DR_gene[rownames(Clinical_tt), ], Clinical_tt$death_censored_graft_loss_event, sd), tapply(DR_nongene[rownames(Clinical_tt), ], Clinical_tt$death_censored_graft_loss_event, sd)))

	        Result_all[nn, ] <- c(Result[1, ], Result[2, ])
	
}

write.table(Result_all, file = paste0(args[1], "_Result.all.xls"), sep = "\t", quote = FALSE, col.names = NA)

