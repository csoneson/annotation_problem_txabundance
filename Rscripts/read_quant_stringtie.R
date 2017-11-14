## Read StringTie gtf

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))

read_quant <- function(file, avefraglength) {
  txdb <- makeTxDbFromGFF(file, format = "gtf")
  ebt <- exonsBy(txdb, "tx")
  gtf <- import(file) 
  gtf <- subset(gtf, type == "transcript")
  mcols(gtf)$length <- width(gtf)
  as.data.frame(mcols(gtf)[, c("transcript_id", "TPM", "length", "cov")]) %>%
    dplyr::rename(transcript = transcript_id) %>%
    dplyr::mutate(TPM = as.numeric(as.character(TPM)),
                  cov = as.numeric(as.character(cov))) %>%
    dplyr::mutate(count = cov * length / avefraglength) %>%
    dplyr::select(transcript, count, TPM)
}