
tableF <- function(tt){

	tmp <- table(tt)[c("0", "1")]

	tmp[is.na(tmp)] <- 0

	tmp

}

biCat <- function(tt){

	Tmp <- table(tt)

	paste0(Tmp["1"], "/", sum(Tmp, na.rm = T), " (", round(Tmp["1"]/(sum(Tmp, na.rm = T))*100, 1), ")")

}

multiCat <- function(tt){

        Tmp <- table(tt)

	Resu <- paste0(Tmp, "/", sum(Tmp, na.rm = T), " (", round(Tmp/sum(Tmp, na.rm = T)*100, 1), ")")

	names(Resu) <- names(Tmp)

	Resu

}

numSum <- function(tt){

	paste0(round(mean(tt, na.rm = T), 2),  "+-", round(sd(tt, na.rm = T), 2))

}

Args <- commandArgs(trailingOnly = T)

Tab <- read.delim(file = "Clinical_Haplo.xls", row.names = 1, as.is = T)

Tab$graft_loss_time_update <- Tab$graft_loss_time_update/30

if(Args[1] == "CC"){

	Tab <- Tab[Tab$R2R == "Caucasian > Caucasian", ]

}else if(Args[1] == "noCC"){

	Tab <- Tab[Tab$R2R != "Caucasian > Caucasian", ]

}


APOL1Geno <- read.delim(file = "../../6.APOL1/GOCAR.xls", row.names = 1, as.is = T)

rownames(APOL1Geno) <- paste0(rownames(APOL1Geno), "Rb")

Tab$APOL1Geno <- APOL1Geno[rownames(Tab), "Rec_APOL1_QC"]

Tab$Rec_APOL1RiskN <- APOL1Geno[rownames(Tab), "Rec_APOL1RiskN"]

Tab$Rec_APOL1RiskN[Tab$Rec_APOL1RiskN == "Unknown"] <- NA


#Tab$GenoGroup2 <- Tab$GenoGroup

#Tab$GenoGroup2[Tab$GenoGroup == "Maj" & Tab$DR_Haplo1_MM == 2] <- "Maj2"
#
#Tab$GenoGroup2[Tab$GenoGroup == "Maj" & Tab$DR_Haplo1_MM == 1] <- "Maj1"
#
#Tab$GenoGroup2[Tab$GenoGroup == "Min" & Tab$DR_Haplo1_MM == 2] <- "Min2"
#
#Tab$GenoGroup2[Tab$GenoGroup == "Min" & Tab$DR_Haplo1_MM == 1] <- "Min1"


Tab$Group <- "Other"

Tab$Group[Tab$GenoGroup == "Maj" & Tab$NEJGenoGroup != "Risk"] <- "Maj"

Tab$Group[Tab$GenoGroup == "Min" & Tab$NEJGenoGroup != "Risk"] <- "Min"

Tab$Group[Tab$NEJGenoGroup == "Risk"] <- "NEJRisk"

Tab$Group <- factor(Tab$Group, levels = c("Other", "NEJRisk", "Maj", "Min"))


#Tab$ExprG[Tab$ExprG == "Medium"] <- "Low"
#
#Tab$transplant_vintage <- Tab$graft_loss_time_update/30
#
#Tab$ExprG <- factor(Tab$ExprG, levels = c("Low", "High"))

Tab$HLA.mm.sum <- factor(Tab$HLA.mm.sum)


Bivariate <- c("death_censored_graft_loss_event", "ACR_event_w_Borderline", "ACR_event_wo_Borderline", "Recurrent_ACR_w_borderline", "Recurrent_ACR_wo_borderline", "Gender", "Donor_Gender", "Induction2", "Donor_Status2")

MutiVariate <- c("Genetic_Rec_Race", "Genetic_Donor_Race", "HLA.mm.sum", "Rec_APOL1RiskN")

Numeric <- c("graft_loss_time_update", "Age", "Donor_Age", "noLIMS1", "Rec_pAFR", "Donor_pAFR")

## Condition

Factor <- Tab$Group

StatBi <- t(apply(Tab[, Bivariate, drop = F], 2, function(x) {x <- factor(x); c(biCat(x), tapply(x, Factor, biCat), round(fisher.test(Reduce(cbind, tapply(x, Factor, tableF)))$p.value, 2))}))

colnames(StatBi) <- c("All", levels(Factor), "Pval")

#StatBiCat <- Reduce(rbind, lapply(BivariateCat, function(x) {x <- factor(x); cbind(multiCat(factor(Tab[[x]])), Reduce(cbind, tapply(factor(Tab[[x]]), Factor, multiCat)), round(fisher.test(Reduce(cbind, tapply(factor(Tab[[x]]), Factor, table)))$p.value, 2))}))

#colnames(StatBi) <- c("All", levels(Factor), "Pval")

StatMulti <- Reduce(rbind, lapply(MutiVariate, function(x) {Tab[[x]] <- factor(Tab[[x]]); cbind(multiCat(Tab[[x]]), Reduce(cbind, tapply(Tab[[x]], Factor, multiCat)), round(chisq.test(Reduce(cbind, tapply(Tab[[x]], Factor, table)))$p.value, 2))}))

colnames(StatMulti) <- c("All", levels(Factor), "Pval")

#StatNum <- t(apply(Tab[, Numeric, drop = F], 2, function(x) {Tmp <- data.frame(Val = x, Fac = Factor); c(numSum(x), tapply(x, Factor, numSum), round(t.test(Val ~ Fac, data = Tmp)$p.value, 2))}))

StatNum <- t(apply(Tab[, Numeric, drop = F], 2, function(x) {x <- as.numeric(x); Tmp <- data.frame(Val = x, Fac = Factor); c(numSum(x), tapply(x, Factor, numSum), round(summary(aov(Val ~ Fac, data = Tmp))[[1]][["Pr(>F)"]][1], 2))}))

colnames(StatNum) <- c("All", levels(Factor), "Pval")


Stat <- rbind(StatBi, StatMulti, StatNum)

write.table(Stat, file = paste0(Args[1], "_Stat_condition.2.xls"), sep = "\t", quote = F, col.names = NA)

