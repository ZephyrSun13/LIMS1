
library(survival)
library(survminer)
library(RColorBrewer)
library(ggpubr)


Colors <- read.delim(file = "ColorLibrary", row.names = 1, as.is = T, header = F)

Clinical <- read.delim(file = "Clinical_Haplo.xls", row.names = 1, as.is = T)

Clinical$GL_time <- Clinical$GL_time/30

Clinical$GenoGroupNEJHaplo <- Clinical$GenoGroup

Clinical$GenoGroupNEJHaplo[which(Clinical$NEJGenoGroup == "Risk")] <- "Risk"

Clinical$GenoGroupNEJHaplo[which(is.na(Clinical$NEJGenoGroup))] <- NA


GOCARHLA <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/db_ZS/GOCAR/HLA.mm.xls", row.names = 1, as.is = T)

Clinical$HLA.3 <- GOCARHLA[rownames(Clinical), "HLA.3"]

Clinical$HLA.4 <- GOCARHLA[rownames(Clinical), "HLA.4"]


Clinical.CC <- Clinical[Clinical$R2R == "Caucasian > Caucasian", ]

ResultAll <- matrix(0, ncol = 8, nrow = 1 , dimnames = list(c("blank"), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N")))

Terms <- c("LIMS1", "NEJRiskHeta", "DR_Haplo1_Heta")

for(tt in Terms){

	fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", tt, " + All + Induction2 + Donor_Status2 + HLA.mm.sum")), data = Clinical)

	ResultAll <- rbind(ResultAll, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n))

}

write.table(ResultAll, file = "HLAOri.xls", quote = F, col.names = NA, sep = "\t")


ResultAll <- matrix(0, ncol = 8, nrow = 1 , dimnames = list(c("blank"), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N")))

Terms <- c("LIMS1", "NEJRiskHeta", "DR_Haplo1_Heta")

for(tt in Terms){
        
        fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", tt, " + All + Induction2 + Donor_Status2 + HLA.mm.sum")), data = Clinical.CC)
        
        ResultAll <- rbind(ResultAll, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n))

}

write.table(ResultAll, file = "HLAOri.CC.xls", quote = F, col.names = NA, sep = "\t")


ResultAll <- matrix(0, ncol = 8, nrow = 1 , dimnames = list(c("blank"), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N")))

Terms <- c("LIMS1", "NEJRiskHeta", "DR_Haplo1_Heta")

for(tt in Terms){
        
        fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", tt, " + All + Induction2 + Donor_Status2 + HLA.4")), data = Clinical)
        
        ResultAll <- rbind(ResultAll, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n))

}

write.table(ResultAll, file = "HLA.4.xls", quote = F, col.names = NA, sep = "\t")


ResultAll <- matrix(0, ncol = 8, nrow = 1 , dimnames = list(c("blank"), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N")))

Terms <- c("LIMS1", "NEJRiskHeta", "DR_Haplo1_Heta")

for(tt in Terms){

        fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", tt, " + All + Induction2 + Donor_Status2 + HLA.4")), data = Clinical.CC)

        ResultAll <- rbind(ResultAll, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n))

}

write.table(ResultAll, file = "HLA.4.CC.xls", quote = F, col.names = NA, sep = "\t")

