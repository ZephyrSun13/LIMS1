
library(survival)
library(survminer)
library(RColorBrewer)
library(ggpubr)


Colors <- read.delim(file = "ColorLibrary", row.names = 1, as.is = T, header = F)

Clinical <- read.delim(file = "../../4.CTOT2/4.SurvivalHaplo/Clinical_Haplo.xls", row.names = 1, as.is = T)

rownames(Clinical) <- paste0("CTOT", rownames(Clinical))


Clinical$GL_time <- Clinical$GL_time/30


CTOTHLA <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/db_ZS/CTOT/HLA.mm.xls", row.names = 1, as.is = T)

#rownames(CTOTHLA) <- CTOTHLA$ID

Clinical$HLA.3 <- CTOTHLA[rownames(Clinical), "HLA.3"]

Clinical$HLA.6 <- CTOTHLA[rownames(Clinical), "HLA.6"]


Clinical.CC <- Clinical[Clinical$R2R == "EUR > EUR", ]

ResultAll <- matrix(0, ncol = 8, nrow = 1 , dimnames = list(c("blank"), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N")))

Terms <- c("LIMS1", "DR_Haplo1_Heta")

for(tt in Terms){

	fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", tt, " + All + induct + DTYPE + HLAmm.6variable")), data = Clinical)

	ResultAll <- rbind(ResultAll, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n))

}

write.table(ResultAll, file = "HLAOri.CTOT.xls", quote = F, col.names = NA, sep = "\t")


ResultAll <- matrix(0, ncol = 8, nrow = 1 , dimnames = list(c("blank"), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N")))

Terms <- c("LIMS1", "DR_Haplo1_Heta")

for(tt in Terms){
        
        fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", tt, " + All + induct + DTYPE + HLAmm.6variable")), data = Clinical.CC)
        
        ResultAll <- rbind(ResultAll, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n))

}

write.table(ResultAll, file = "HLAOri.CTOT.CC.xls", quote = F, col.names = NA, sep = "\t")


ResultAll <- matrix(0, ncol = 8, nrow = 1 , dimnames = list(c("blank"), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N")))


for(tt in Terms){
        
        fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", tt, " + All + induct + DTYPE + HLA.6")), data = Clinical)
        
        ResultAll <- rbind(ResultAll, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n))

}

write.table(ResultAll, file = "HLA.6.CTOT.xls", quote = F, col.names = NA, sep = "\t")


ResultAll <- matrix(0, ncol = 8, nrow = 1 , dimnames = list(c("blank"), c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N")))


for(tt in Terms){

        fit <- coxph(as.formula(paste0("Surv(GL_time, GL_event) ~", tt, " + All + induct + DTYPE + HLA.6")), data = Clinical.CC)

        ResultAll <- rbind(ResultAll, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n))

}

write.table(ResultAll, file = "HLA.6.CTOT.CC.xls", quote = F, col.names = NA, sep = "\t")

