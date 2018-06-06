################################################################################
##                                                                            ##
## Help function to read abundance.tsv file from kallisto                     ##
##                                                                            ##
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

read_quant <- function(file, avefraglength) {
  tmp <- read.delim(file, header = TRUE, as.is = TRUE)
  idx <- grep("^STRG\\.|^CHS\\.", tmp$target_id, invert = TRUE)
  tmp$target_id[idx] <- gsub("\\.[0-9]+$", "", tmp$target_id[idx])
  tmp %>% dplyr::rename(transcript = target_id, count = est_counts, TPM = tpm) %>%
    dplyr::select(transcript, count, TPM)
}