
Args <- commandArgs(trailingOnly = T)

library(survival)
library(survminer)
library(RColorBrewer)
library(ggpubr)


Colors <- read.delim(file = "ColorLibrary", row.names = 1, as.is = T, header = F)

Clinical <- read.delim(file = "../../4.CTOT2/4.SurvivalHaplo/Clinical_Haplo.xls", row.names = 1, as.is = T)

Clinical$GL_time <- Clinical$GL_time/30

Clinical2 <- read.delim(file = "../../6.APOL1/CTOT01.worksheet.txt", sep = "\t", row.names = 2, stringsAsFactors = FALSE)

Clinical$ACR_borderline <- Clinical2[rownames(Clinical), "All.AR"]

Clinical$ACR_without_borderline <- Clinical2[rownames(Clinical), "AR.1A"]


#Clinical$GenoGroupNEJHaplo <- Clinical$GenoGroup

#Clinical$GenoGroupNEJHaplo[which(Clinical$NEJGenoGroup == "Risk")] <- "Risk"

#Clinical$GenoGroupNEJHaplo[which(is.na(Clinical$NEJGenoGroup))] <- NA


Clinical$AR <- Clinical[[Args[1]]]

Clinical.CC <- Clinical[Clinical$R2R == "EUR > EUR", ]

ResultAll <- matrix(0, ncol = 8, nrow = 1 , dimnames = list(c("blank"), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N")))

Terms <- c("LIMS1", "DR_Haplo1_Heta")

for(tt in Terms){

	fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", tt, " + All + induct + DTYPE + HLAmm.6variable + AR")), data = Clinical)

	ResultAll <- rbind(ResultAll, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n))

}

write.table(ResultAll, file = paste0(Args[1], ".CTOT.corr.xls"), quote = F, col.names = NA, sep = "\t")


ResultAll <- matrix(0, ncol = 8, nrow = 1 , dimnames = list(c("blank"), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N")))

Terms <- c("LIMS1", "DR_Haplo1_Heta")

for(tt in Terms){
        
        fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", tt, " + All + induct + DTYPE + HLAmm.6variable + AR")), data = Clinical.CC)
        
        ResultAll <- rbind(ResultAll, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n))

}

write.table(ResultAll, paste0(file = Args[1], ".CC.CTOT.corr.xls"), quote = F, col.names = NA, sep = "\t")

