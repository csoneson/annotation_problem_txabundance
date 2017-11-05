args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(inrds)
print(outtxt)

x <- readRDS(inrds)
write.table(x$gene, file = outtxt, row.names = FALSE, col.names = FALSE, 
            sep = "\t", quote = FALSE)

sessionInfo()
date()
