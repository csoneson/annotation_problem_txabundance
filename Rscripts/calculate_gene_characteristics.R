args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(txome)
print(outrds)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(dplyr)
  library(ggplot2)
})

gtf <- rtracklayer::import(gtf, format = "gtf")
txome <- readDNAStringSet(txome)

## Number of exons per transcript, total exon length per transcript
exons <- as.data.frame(gtf) %>% dplyr::filter(type == "exon") %>% 
  dplyr::group_by(transcript_id) %>%
  dplyr::summarize(nbr_exons = length(exon_id),
                   mean_exon_length = mean(width),
                   longest_exon = max(width),
                   shortest_exon = min(width),
                   median_exon_length = median(width),
                   tx_length = sum(width),
                   gene_id = gene_id[1])

## 3' UTRs per transcript
utrs <- as.data.frame(gtf) %>% dplyr::filter(type == "three_prime_utr") %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarize(start = min(start), end = max(end), width = sum(width), 
                   seqnames = seqnames[1], gene_id = gene_id[1], strand = strand[1])

## If a gene has multiple 3'UTRs starting in the same place but with different
## length, get the length difference
tmp <- utrs %>% dplyr::select(-transcript_id)
tmp <- split(tmp, tmp$gene_id)
tmp2 <- sapply(tmp, function(w) {
  if (all(w$strand == "-")) {
    w %>% group_by(end) %>% summarize(widthdiff = max(width) - min(width)) %>% 
      select(widthdiff) %>% max
  } else if (all(w$strand == "+")) {
    w %>% group_by(start) %>% summarize(widthdiff = max(width) - min(width)) %>% 
      select(widthdiff) %>% max
  } else {
    0
  }
})

## Summarize on gene level
gene_char <- dplyr::full_join(
  exons %>% dplyr::group_by(gene_id) %>%
    dplyr::summarize(nbr_transcripts_gtf = length(transcript_id),
                     ave_transcript_length_gtf = mean(tx_length),
                     min_transcript_length_gtf = min(tx_length),
                     max_transcript_length_gtf = max(tx_length),
                     median_transcript_length_gtf = median(tx_length),
                     ave_nbr_exons_per_tx = mean(nbr_exons),
                     min_nbr_exons_per_tx = min(nbr_exons),
                     max_nbr_exons_per_tx = max(nbr_exons),
                     median_nbr_exons_per_tx = median(nbr_exons)),
  utrs %>% dplyr::select(-transcript_id) %>% dplyr::distinct() %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::summarize(nbr_distinct_3putrs = length(unique(paste0(start, ".", end))),
                     max_3putr_length = max(width),
                     min_3putr_length = min(width),
                     ave_3putr_length = mean(width),
                     median_3putr_length = median(width))
) %>% dplyr::full_join(data.frame(gene_id = names(tmp2), 
                                  length_diff_3putrs_samestart = tmp2,
                                  stringsAsFactors = FALSE))

## Add information about the genes that are not present in the gtf
tx_info <- data.frame(width = width(txome), id = names(txome),
                      stringsAsFactors = FALSE) %>%
  dplyr::mutate(transcript_id = sapply(strsplit(id, " "), .subset, 1)) %>%
  dplyr::mutate(gene_id = sapply(strsplit(id, " "), 
                                 function(w) gsub("^gene:", "", w[grep("^gene:", w)]))) %>%
  dplyr::select(transcript_id, gene_id, width) %>%
  dplyr::mutate(transcript_id = gsub("\\.[0-9]+$", "", transcript_id),
                gene_id = gsub("\\.[0-9]+$", "", gene_id))
gene_info <- tx_info %>% dplyr::group_by(gene_id) %>%
  dplyr::summarize(nbr_transcripts_fasta = length(transcript_id),
                   ave_transcript_length_fasta = mean(width),
                   min_transcript_length_fasta = min(width),
                   max_transcript_length_fasta = max(width),
                   median_transcript_length_fasta = median(width))

## Whenever a gene is present in both fasta and gtf, the information agrees.
## Generate new columns consolidating the information from both sources.
gene_info <- dplyr::full_join(gene_char, gene_info) %>%
  dplyr::mutate(nbr_transcripts = pmax(nbr_transcripts_gtf, nbr_transcripts_fasta, na.rm = TRUE),
                ave_transcript_length = pmax(ave_transcript_length_gtf,
                                             ave_transcript_length_fasta, na.rm = TRUE),
                min_transcript_length = pmax(min_transcript_length_gtf,
                                             min_transcript_length_fasta, na.rm = TRUE),
                max_transcript_length = pmax(max_transcript_length_gtf,
                                             max_transcript_length_fasta, na.rm = TRUE),
                median_transcript_length = pmax(median_transcript_length_gtf,
                                                median_transcript_length_fasta, na.rm = TRUE))


saveRDS(gene_info, file = outrds)
sessionInfo()
date()
