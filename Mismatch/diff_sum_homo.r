
args <- commandArgs(trailingOnly = TRUE)

#DR <- read.table(file = args[1], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

DR <- read.table(file = args[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

DR[is.na(DR)] <- 0

DR[DR==1] <- 0

DR[DR==2] <- 1

#Ann <- read.table(file = args[2], )

DR_freq <- colSums(DR)

Mut_freq <- rowSums(DR)/ncol(DR)

#DR_freq <- DR_freq[order(DR_freq, decreasing = TRUE)]

write.table(DR_freq, file = paste0(args[1], ".homo.freq.xls"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

write.table(DR[Mut_freq >= 0.05, ], file = paste0(args[1], ".comm.homo"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

write.table(colSums(DR[Mut_freq >= 0.05, ]), file = paste0(args[1], ".comm.homo.freq.xls"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

write.table(colSums(DR[Mut_freq < 0.05, ]), file = paste0(args[1], ".rare.homo.freq.xls"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

write.table(Mut_freq, file = paste0(args[1], ".homo.mutfreq.xls"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

#Stat <- data.frame(Sample = names(DR_freq), freq = DR_freq)

#library(ggplot2)
#
#p <- ggplot(Stat, aes(x=freq)) +
#        geom_density() +
#        geom_vline(aes(xintercept=mean(freq)), color="blue", linetype="dashed", size=1) +
#        theme_minimal()
#
#ggsave(p, file = paste0(args[1], ".Frequency.pdf"))

