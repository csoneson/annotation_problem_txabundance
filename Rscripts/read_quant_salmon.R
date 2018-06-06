################################################################################
##                                                                            ##
## Help function to read quant.sf file from Salmon                            ##
##                                                                            ##
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

read_quant <- function(file, avefraglength) {
  tmp <- read.delim(file, header = TRUE, as.is = TRUE)
  idx <- grep("^STRG\\.|^CHS\\.", tmp$Name, invert = TRUE)
  tmp$Name[idx] <- gsub("\\.[0-9]+$", "", tmp$Name[idx])
  tmp %>% dplyr::rename(transcript = Name, count = NumReads) %>%
    dplyr::select(transcript, count, TPM)
}