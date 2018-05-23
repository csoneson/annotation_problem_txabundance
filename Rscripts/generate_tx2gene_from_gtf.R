################################################################################
##                                                                            ##
## Generate transcript-to-gene conversion table from gtf file                 ##
##                                                                            ##
## Inputs:                                                                    ##
## * gtf: gtf file                                                            ##
## * outrds: output file                                                      ##
##                                                                            ##
## Outputs:                                                                   ##
## * A transcript-to-gene conversion table                                    ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(outrds)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(dplyr)
})

gtf <- import(gtf)

gtfsub <- as.data.frame(gtf) %>% dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::rename(tx = transcript_id, 
                gene = gene_id)


saveRDS(gtfsub, file = outrds)

sessionInfo()
date()
