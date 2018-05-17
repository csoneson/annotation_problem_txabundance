args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(txfasta)
print(outfasta)
print(readdir)
print(readbasename)
print(readlen)

suppressPackageStartupMessages({
  library(Biostrings)
  library(rtracklayer)
  library(dplyr)
  library(polyester)
  library(ShortRead)
})

gtf <- rtracklayer::import(gtf, format = "gtf")
txfasta <- Biostrings::readDNAStringSet(txfasta)
names(txfasta) <- gsub("\\.[0-9]+$", "", sapply(strsplit(names(txfasta), " "), .subset, 1))

## Get the coordinates of the 3'UTR for each transcript
utrs <- as.data.frame(gtf) %>% dplyr::filter(type == "three_prime_utr") %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarize(start = min(start), end = max(end), width = sum(width), 
                   seqnames = seqnames[1], gene_id = gene_id[1], strand = strand[1])

## If a gene has multiple 3'UTRs starting in the same place but with different
## length, get the length difference
tmp <- utrs %>% dplyr::select(-transcript_id)
tmp <- split(tmp, tmp$gene_id)
utr_length_diff <- sapply(tmp, function(w) {
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

## Get the names of the genes with large differences in 3'UTR lengths
genes_of_interest <- names(utr_length_diff)[utr_length_diff > 1000]
length(genes_of_interest)

## Initialize a DNAStringSet to hold the modified transcripts
modified_transcripts <- DNAStringSet()

## For each of these genes, find the 3'UTRs that start in the same place but
## have different lengths, and get the transcripts corresponding to the shortest
## and longest of these 3'UTRs. Then generate a new sequence representing one of
## these transcripts, but with the 3'UTR replaced with that of the other
## transcript.
set.seed(42)
for (g in genes_of_interest) {
  ## Get all 3'UTRs for the gene
  utrsub <- subset(utrs, gene_id == g)
  
  ## Group 3'UTRs by their starting position
  if (all(utrsub$strand == "-")) {
    utrsub <- utrsub %>% dplyr::group_by(end)
  } else {
    utrsub <- utrsub %>% dplyr::group_by(start)
  }
  
  ## Keep only starting positions with at least two 3'UTRs, and among those,
  ## keep only the starting position with the largest difference between the
  ## longest and shortest 3'UTR. For this position, keep the shortest and
  ## longest 3'UTR.
  utrsub <- utrsub %>% dplyr::mutate(n = length(width)) %>% 
    dplyr::filter(n > 1) %>% dplyr::mutate(lengthdiff = max(width) - min(width)) %>%
    dplyr::ungroup() %>% dplyr::filter(lengthdiff == max(lengthdiff)) %>%
    dplyr::arrange(width) %>% dplyr::filter(row_number() %in% c(1, n()))
  stopifnot(nrow(utrsub) == 2)
  
  ## Randomly choose one of the transcripts to retain as the basis, and set the
  ## other one as the transcript to get the 3'UTR from.
  base_tx <- sample(utrsub$transcript_id, 1)
  utr_tx <- setdiff(utrsub$transcript_id, base_tx)
  
  ## Get the transcript sequence of the base transcript and the 3'UTR-providing
  ## transcript
  base_tx_seq <- as.character(txfasta[base_tx])
  utr_tx_seq <- as.character(txfasta[utr_tx])
  
  ## Get the 3'UTR lengths for the two transcripts
  base_tx_utr_width <- utrsub$width[utrsub$transcript_id == base_tx]
  utr_tx_utr_width <- utrsub$width[utrsub$transcript_id == utr_tx]
  
  ## Get the base sequence up until the original 3'UTR
  base_seq <- substring(base_tx_seq, first = 1, last = nchar(base_tx_seq) - base_tx_utr_width)
  
  ## Get the 3'UTR sequence to add
  utr_seq <- substring(utr_tx_seq, first = nchar(utr_tx_seq) - utr_tx_utr_width + 1, last = nchar(utr_tx_seq))
  
  ## Generate the final modified sequence and add to the collection
  modified_transcripts <- c(modified_transcripts,
                            DNAStringSet(x = structure(paste0(base_seq, utr_seq), 
                                                       names = paste0(base_tx, "_utrfrom_", utr_tx))))
}

## Select additional transcripts randomly (excluding genes considered above) and
## assign them counts
available_transcripts <- setdiff(names(txfasta), 
                                 gtf$transcript_id[gtf$gene_id %in% genes_of_interest])
set.seed(42)
additional_transcripts <- sample(available_transcripts, size = 10000, replace = FALSE)
additional_counts <- c(rmultinom(n = 1, size = 10e6, prob = runif(10000)))
## Check the number of isoforms of the genes of the selected transcripts
table(as.data.frame(gtf) %>% dplyr::filter(type == "transcript") %>%
        dplyr::group_by(gene_id) %>% dplyr::mutate(nbr_tx = length(transcript_id)) %>%
        dplyr::filter(transcript_id %in% additional_transcripts) %>%
        dplyr::select(gene_id, nbr_tx) %>% dplyr::distinct() %>% dplyr::pull(nbr_tx))

## Generate final set of transcripts to simulate from and number of reads to generate
transcripts_to_simulate_from <- c(modified_transcripts, txfasta[additional_transcripts])
reads_per_transcript <- c(rep(1000, length(modified_transcripts)),
                          additional_counts)

## Subset to transcripts that are at least as long as the fragment length
idx <- which(width(transcripts_to_simulate_from) >= readlen)
transcripts_to_simulate_from <- transcripts_to_simulate_from[idx]
reads_per_transcript <- reads_per_transcript[idx]

## Write modified and additional transcripts to fasta file
writeXStringSet(transcripts_to_simulate_from, 
                filepath = outfasta)

## Simulate reads with polyester. Generates files sample_01_1.fasta.gz and
## sample_01_2.fasta.gz in the readdir
polyester::simulate_experiment(fasta = outfasta, 
                               outdir = readdir, 
                               fold_changes = 1,
                               num_reps = c(1, 1),
                               reads_per_transcript = reads_per_transcript,
                               size = 500,
                               paired = TRUE,
                               reportCoverage = FALSE,
                               readlen = readlen,
                               distr = "normal",
                               fraglen = 300,
                               fragsd = 25,
                               strand_specific = TRUE,
                               seed = 42,
                               gzip = TRUE)

## Write reads to fastq files and save list of modified genes
saveRDS(genes_of_interest, file = paste0(readdir, "/", readbasename,
                                         "_modified_genes.rds"))

date()
sessionInfo()
