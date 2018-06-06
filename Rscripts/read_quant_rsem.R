################################################################################
##                                                                            ##
## Help function to read isoforms.results file from RSEM                      ##
##                                                                            ##
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

read_quant <- function(file, avefraglength) {
  tmp <- read.delim(file, header = TRUE, as.is = TRUE)
  idx <- grep("^STRG\\.|^CHS\\.", tmp$transcript_id, invert = TRUE)
  tmp$transcript_id[idx] <- gsub("\\.[0-9]+$", "", tmp$transcript_id[idx])
  tmp %>% dplyr::rename(transcript = transcript_id, count = expected_count) %>%
    dplyr::select(transcript, count, TPM)
}