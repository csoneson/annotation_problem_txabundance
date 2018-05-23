################################################################################
##                                                                            ##
## Help function to read bam_count_reads.tsv file from Wub                    ##
##                                                                            ##
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

read_quant <- function(file, avefraglength) {
  read.delim(file, header = TRUE, as.is = TRUE) %>%
    dplyr::select(Reference, Count) %>%
    dplyr::mutate(TPM = Count/sum(Count) * 1e6) %>% ## assume one read/molecule
    dplyr::mutate(Reference = gsub("\\.[0-9]+", "", Reference)) %>%
    dplyr::rename(transcript = Reference, count = Count) %>%
    dplyr::select(transcript, count, TPM)
}
