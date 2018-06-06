################################################################################
##                                                                            ##
## Help function to read bam_count_reads.tsv file from Wub                    ##
##                                                                            ##
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

read_quant <- function(file, avefraglength) {
  tmp <- read.delim(file, header = TRUE, as.is = TRUE) %>%
    dplyr::select(Reference, Count) %>%
    dplyr::mutate(TPM = Count/sum(Count) * 1e6) ## assume one read/molecule
  idx <- grep("^STRG\\.|^CHS\\.", tmp$Reference, invert = TRUE)
  tmp$Reference[idx] <- gsub("\\.[0-9]+$", "", tmp$Reference[idx])
  tmp %>% dplyr::rename(transcript = Reference, count = Count) %>%
    dplyr::select(transcript, count, TPM)
}
