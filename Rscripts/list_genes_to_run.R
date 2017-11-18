args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(inrds)
print(outtxt)

suppressPackageStartupMessages(library(dplyr))

x <- readRDS(inrds) %>% dplyr::filter(maxnbrex > 1 & 
                                        count > 1000 & 
                                        maxlength > 350)
write.table(x$gene, file = outtxt, row.names = FALSE, col.names = FALSE, 
            sep = "\t", quote = FALSE)

sessionInfo()
date()
