args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(inrds)  ## output from alpine_get_predicted_coverage.R
print(tx2gene)
print(outtxt)

suppressPackageStartupMessages(library(dplyr))

x <- readRDS(inrds)
## Keep only transcripts with at least one junction
x <- x[sapply(x, function(w) length(w$junctions)) > 0]
tx2gene <- readRDS(tx2gene) %>% dplyr::mutate(tx = gsub("\\.[0-9]+$", "", tx),
                                              gene = gsub("\\.[0-9]+$", "", gene))
genes <- unique(tx2gene$gene[match(names(x), tx2gene$tx)])

write.table(genes, file = outtxt, row.names = FALSE, col.names = FALSE, 
            sep = "\t", quote = FALSE)

sessionInfo()
date()
