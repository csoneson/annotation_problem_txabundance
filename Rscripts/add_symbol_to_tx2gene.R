################################################################################
##                                                                            ##
## Add the gene symbol to a tx2gene data frame                                ##
##                                                                            ##
## Inputs:                                                                    ##
## * tx2gene: tx2gene data frame                                              ##
## * info: data frame with gene and symbol columns                            ##
## * outrds: output file where the extended tx2gene data frame will be saved  ##
##                                                                            ##
## Outputs:                                                                   ##
## * An extended tx2gene data frame                                           ##
##                                                                            ##
################################################################################

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

info <- info %>% dplyr::select(gene, symbol) %>% dplyr::distinct()
idx <- grep("^STRG\\.|^CHS\\.", info$gene, invert = TRUE)
info$gene[idx] <- gsub("\\.[0-9]+$", "", info$gene[idx])

tx2gene <- tx2gene %>% 
  dplyr::left_join(info, by = "gene") %>%
  dplyr::mutate(symbol = replace(symbol, is.na(symbol), gene[is.na(symbol)]))

saveRDS(tx2gene, outrds)

date()
sessionInfo()
