## Read StringTie gtf

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(dplyr))

read_quant <- function(file, avefraglength) {
  gtf <- import(file)
  gtfex <- subset(gtf, type == "exon")
  txlengths <- as.data.frame(gtfex) %>% dplyr::group_by(transcript_id) %>%
    dplyr::summarize(width = sum(width))
  gtftx <- subset(gtf, type == "transcript")
  as.data.frame(gtftx) %>% dplyr::select(transcript_id, TPM, cov) %>%
    dplyr::left_join(txlengths) %>%
    dplyr::rename(transcript = transcript_id) %>%
    dplyr::mutate(TPM = as.numeric(as.character(TPM)),
                  cov = as.numeric(as.character(cov))) %>%
    dplyr::mutate(count = cov * width / avefraglength) %>%
    dplyr::select(transcript, count, TPM)
}