## Read kallisto abundance.tsv

read_quant <- function(file, avefraglength) {
  read.delim(file, header = TRUE, as.is = TRUE) %>%
    dplyr::mutate(target_id = gsub("\\.[0-9]+", "", target_id)) %>%
    dplyr::rename(transcript = target_id, count = est_counts, TPM = tpm) %>%
    dplyr::select(transcript, count, TPM)
}