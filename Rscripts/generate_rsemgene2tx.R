args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(tx2gene)
print(rsemgene2tx)

x <- readRDS(tx2gene)

write.table(x[, c("gene", "tx")], file = rsemgene2tx, 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

sessionInfo()
date()
