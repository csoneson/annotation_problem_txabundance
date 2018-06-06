################################################################################
##                                                                            ##
## Calculate gene characteristics for CHESS annotation                        ##
##                                                                            ##
## Inputs:                                                                    ##
## * gtf: gtf file                                                            ##
## * outrds: output file                                                      ##
##                                                                            ##
## Outputs:                                                                   ##
## * A data frame with gene characteristics                                   ##
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
  library(Biostrings)
  library(dplyr)
  library(ggplot2)
})

gtf <- rtracklayer::import(gtf, format = "gtf")

## Number of exons per transcript, total exon length per transcript
exons <- as.data.frame(gtf) %>% dplyr::filter(type == "exon") %>% 
  dplyr::group_by(transcript_id) %>%
  dplyr::summarize(nbr_exons = length(exon_id),
                   mean_exon_length = mean(width),
                   longest_exon = max(width),
                   shortest_exon = min(width),
                   median_exon_length = median(width),
                   tx_length = sum(width),
                   gene_id = gene_id[1])

## Summarize on gene level
gene_info <- exons %>% dplyr::group_by(gene_id) %>%
    dplyr::summarize(nbr_transcripts = length(transcript_id),
                     ave_transcript_length = mean(tx_length),
                     min_transcript_length = min(tx_length),
                     max_transcript_length = max(tx_length),
                     median_transcript_length = median(tx_length),
                     ave_nbr_exons_per_tx = mean(nbr_exons),
                     min_nbr_exons_per_tx = min(nbr_exons),
                     max_nbr_exons_per_tx = max(nbr_exons),
                     median_nbr_exons_per_tx = median(nbr_exons))


saveRDS(gene_info, file = outrds)
sessionInfo()
date()
