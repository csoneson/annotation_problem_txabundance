################################################################################
##                                                                            ##
## Help function to read isoforms.results file from RSEM                      ##
##                                                                            ##
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

read_quant <- function(file, avefraglength) {
  read.delim(file, header = TRUE, as.is = TRUE) %>%
    dplyr::mutate(transcript_id = gsub("\\.[0-9]+", "", transcript_id)) %>%
    dplyr::rename(transcript = transcript_id, count = expected_count) %>%
    dplyr::select(transcript, count, TPM)
}