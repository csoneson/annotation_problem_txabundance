args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(tx2gene)
print(info)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
})

tx2gene <- readRDS(tx2gene)
info <- readRDS(info)

tx2gene <- tx2gene %>% 
  dplyr::left_join(info %>% dplyr::select(gene, symbol) %>%
                     dplyr::mutate(gene = gsub("\\.[0-9]+$", "", gene)), by = "gene") %>%
  dplyr::mutate(symbol = replace(symbol, is.na(symbol), gene[is.na(symbol)]))

saveRDS(tx2gene, outrds)

date()
sessionInfo()
