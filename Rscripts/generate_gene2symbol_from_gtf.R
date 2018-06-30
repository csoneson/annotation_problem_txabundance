################################################################################
##                                                                            ##
## Create mapping between StringTie ID and Ensembl ID                         ##
##                                                                            ##
## Inputs:                                                                    ##
## * gtf: StringTie gtf file                                                  ##
## * outrds: output file where the gene2symbol data frame will be saved       ##
##                                                                            ##
## Outputs:                                                                   ##
## * A data frame mapping StringTie IDs to Ensembl IDs                        ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(rtracklayer)
})

gtf <- import(gtf)
gtf <- subset(gtf, type == "transcript")
gene2symbol <- data.frame(gene = gtf$gene_id, 
                          symbol = gtf$ref_gene_id,
                          stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(symbol)) %>% 
  dplyr::distinct() %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(symbol = paste(symbol, collapse = "__")) %>%
  as.data.frame()
gene2symbol <- rbind(gene2symbol, 
                     data.frame(gene = setdiff(gtf$gene_id, gene2symbol$gene),
                                symbol = NA,
                                stringsAsFactors = FALSE))

saveRDS(gene2symbol, file = outrds)

date()
sessionInfo()