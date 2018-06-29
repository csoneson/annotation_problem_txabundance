################################################################################
##                                                                            ##
## Create mapping between StringTie ID and Ensembl ID for transcripts         ##
##                                                                            ##
## Inputs:                                                                    ##
## * gtf: StringTie gtf file                                                  ##
## * outrds: output file where the tx2symbol data frame will be saved         ##
##                                                                            ##
## Outputs:                                                                   ##
## * A data frame mapping StringTie transcript IDs to Ensembl IDs             ##
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
tx2symbol <- data.frame(tx = gtf$transcript_id, 
                        symbol = gtf$reference_id,
                        stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(symbol)) %>% 
  dplyr::distinct()
tx2symbol <- rbind(tx2symbol, 
                   data.frame(tx = setdiff(gtf$transcript_id, tx2symbol$tx),
                              symbol = NA,
                              stringsAsFactors = FALSE))

saveRDS(tx2symbol, file = outrds)

date()
sessionInfo()