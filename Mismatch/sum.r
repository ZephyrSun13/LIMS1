
Args <- commandArgs(trailingOnly = TRUE)

Tab <- read.table(file = Args[1], header = FALSE, sep = "\t")

Group <- Tab[[1]]

Tab <- Tab[, -1]

SS <- apply(Tab, 2, function(x) tapply(x, Group, sum))

write.table(SS, file = paste0(Args[1], ".Gene"), sep = "\t", quote = FALSE, col.names = FALSE)

