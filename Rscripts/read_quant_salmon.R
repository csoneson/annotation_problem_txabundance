## Read Salmon quant.sf file

suppressPackageStartupMessages(library(dplyr))

read_quant <- function(file, avefraglength) {
  read.delim(file, header = TRUE, as.is = TRUE) %>%
    dplyr::mutate(Name = gsub("\\.[0-9]+", "", Name)) %>%
    dplyr::rename(transcript = Name, count = NumReads) %>%
    dplyr::select(transcript, count, TPM)
}