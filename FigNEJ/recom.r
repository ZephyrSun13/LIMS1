
Tab <- read.delim(file = "Haplo.heat2.1000G.xls", row.names = 1, as.is = T)

List <- list()

for(nn in 1:(ncol(Tab)/2)){

	List[[nn]] <- table(factor(Tab[, nn*2] - Tab[, nn*2-1], levels = c(-1, 0, 1)))

}

Sum <- Reduce(rbind, List)

Recom <- which(Sum[, 1] >= 3 & Sum[, 3] >= 3)

length(Recom)

Tab[, sort(c(Recom * 2 - 1, Recom * 2))]

write.table(Tab[, sort(c(Recom * 2 - 1, Recom * 2))], file = "1000GRecom.xls", sep = "\t", quote = F, col.names = F)

Sams <- unlist(read.table(file = "Sam2504", as.is = T))

Race <- read.delim(file = "integrated_call_samples_v3.20130502.ALL.panel", row.names = 1, as.is = T)

table(Race$super_pop)

table(Race[Sams[Recom], "super_pop"])

table(Race[Sams[Recom], "super_pop"])/table(Race$super_pop)[names(table(Race[Sams[Recom], "super_pop"]))]


Tab <- read.delim(file = "Haplo.heat2.GoCAR.xls", row.names = 1, as.is = T)

List <- list()

for(nn in 1:(ncol(Tab)/2)){

        List[[nn]] <- table(factor(Tab[, nn*2] - Tab[, nn*2-1], levels = c(-1, 0, 1)))

}

Sum <- Reduce(rbind, List)

Recom <- which(Sum[, 1] >= 3 & Sum[, 3] >= 3)

length(Recom)


Tab <- read.delim(file = "Haplo.heat2.CTOT.xls", row.names = 1, as.is = T)

List <- list()

for(nn in 1:(ncol(Tab)/2)){

        List[[nn]] <- table(factor(Tab[, nn*2] - Tab[, nn*2-1], levels = c(-1, 0, 1)))

}

Sum <- Reduce(rbind, List)

Recom <- which(Sum[, 1] >= 3 & Sum[, 3] >= 3)

length(Recom)

