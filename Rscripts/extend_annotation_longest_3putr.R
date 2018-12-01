################################################################################
##                                                                            ##
## Extend short 3'UTRs to the longest 3'UTR starting in the same position     ##
##                                                                            ##
## Inputs:                                                                    ##
## * gtffile: gtf file                                                        ##
## * txfastafile: transcriptome fasta file                                    ##
## * outbase: base name of output files (gtf, fasta, tx2gene)                 ##
##                                                                            ##
## Outputs:                                                                   ##
## * Expanded gtf and fasta files, plus tx2gene file                          ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtffile)
print(txfastafile)
print(outbase)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(dplyr)
  library(GenomicFeatures)
  library(BSgenome)
  library(BSgenome.Hsapiens.NCBI.GRCh38)
})

## Read reference files
gtf <- rtracklayer::import(gtffile, format = "gtf")
txdb <- makeTxDbFromGFF(gtffile)
ebt <- exonsBy(txdb, by = "tx", use.names = TRUE)
txfasta <- Biostrings::readDNAStringSet(txfastafile)
names(txfasta) <- gsub("\\.[0-9]+$", "", sapply(strsplit(names(txfasta), " "), .subset, 1))

## ========================================================================== ##
## Figure out which transcripts should get a new 3'UTR (and which 3'UTR that should be)
## ========================================================================== ##
## Get the coordinates of the 3'UTR for each transcript
utrs <- subset(gtf, type == "three_prime_utr")

## Get the total length of each 3'UTR
utrdf <- as.data.frame(utrs) %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::summarize(start = min(start), end = max(end), width = sum(width), 
                   seqnames = seqnames[1], gene_id = gene_id[1], 
                   strand = strand[1], gene_version = gene_version[1],
                   gene_name = gene_name[1], gene_source = gene_source[1],
                   gene_biotype = gene_biotype[1])

## If a gene has multiple 3'UTRs starting in the same place but with different
## length, get the length difference
tmp <- utrdf %>% dplyr::select(-transcript_id)
tmp <- split(tmp, tmp$gene_id)
utr_length_diff <- sapply(tmp, function(w) {
  if (all(w$strand == "-")) {
    w %>% dplyr::group_by(end) %>% 
      dplyr::summarize(widthdiff = max(width) - min(width)) %>% 
      dplyr::select(widthdiff) %>% max
  } else if (all(w$strand == "+")) {
    w %>% dplyr::group_by(start) %>% 
      dplyr::summarize(widthdiff = max(width) - min(width)) %>% 
      dplyr::select(widthdiff) %>% max
  } else {
    0
  }
})

## Get the names of the genes with at least one 3'UTR start position with 3'UTRs
## of different lengths
genes_of_interest <- names(utr_length_diff)[utr_length_diff > 0]
length(genes_of_interest)

## For each of these genes, find the 3'UTRs that start in the same place but
## have different lengths.
utr_exchange <- do.call(dplyr::bind_rows, lapply(genes_of_interest, function(g) {
  ## Get all 3'UTRs for the gene
  utrsub <- subset(utrdf, gene_id == g)
  
  ## Group 3'UTRs by their starting position. If there are multiple UTRs with
  ## the same starting position and equal (=longest) length, keep only one, to
  ## avoid ambiguities when replacing the other transcripts' UTRs, and to avoid
  ## creating identical transcripts.
  if (all(utrsub$strand == "-")) {
    utrsub <- utrsub %>% dplyr::group_by(end) %>%
      dplyr::mutate(maxwidth = max(width)) %>%
      dplyr::group_by(end, maxwidth) %>% 
      dplyr::mutate(idx = seq_len(length(start))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(width != maxwidth | idx == 1) %>%
      dplyr::group_by(end)
  } else {
    utrsub <- utrsub %>% dplyr::group_by(start) %>%
      dplyr::mutate(maxwidth = max(width)) %>%
      dplyr::group_by(start, maxwidth) %>% 
      dplyr::mutate(idx = seq_len(length(end))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(width != maxwidth | idx == 1) %>%
      dplyr::group_by(start)
  }

  ## Keep only starting positions with at least two 3'UTRs, and generate the 
  ## summary table of which transcripts should get a new 3'UTR
  utrsub %>% dplyr::mutate(n = length(width)) %>% 
    dplyr::filter(n > 1) %>% 
    dplyr::mutate(donating_tx = transcript_id[which.max(width)]) %>%
    dplyr::rename(receiving_tx = transcript_id) %>%
    dplyr::filter(receiving_tx != donating_tx) %>%
    dplyr::ungroup() %>%
    dplyr::select(receiving_tx, donating_tx, gene_id, gene_version,
                  gene_name, gene_source, gene_biotype)
}))

## ========================================================================== ##
## Expand ebt with new transcripts
## ========================================================================== ##
idx <- seq_len(nrow(utr_exchange))
names(idx) <- paste0(utr_exchange$receiving_tx, "longUTR")
ebt_expanded <- as(lapply(idx, function(i) {
  tmp1 <- ebt[[utr_exchange$receiving_tx[i]]]
  mcols(tmp1) <- NULL
  tmp2 <- subset(utrs, transcript_id == utr_exchange$donating_tx[i])
  mcols(tmp2) <- NULL
  tmp <- GenomicRanges::reduce(c(tmp1, tmp2))
  ## Important that the exons are in order for the sequence extraction later
  if (all(strand(tmp) == "-")) {
    tmp <- tmp[order(start(tmp), decreasing = TRUE)]
  } else {
    tmp <- tmp[order(start(tmp), decreasing = FALSE)]
  }
  mcols(tmp) <- data.frame(
    source = "local",
    type = "exon",
    score = NA_real_,
    phase = NA_integer_,
    gene_id = utr_exchange$gene_id[i],
    gene_version = utr_exchange$gene_version[i],
    gene_name = utr_exchange$gene_name[i],
    gene_source = utr_exchange$gene_source[i],
    gene_biotype = utr_exchange$gene_biotype[i],
    transcript_id = paste0(utr_exchange$receiving_tx[i], "longUTR"),
    transcript_version = 1,
    transcript_name = paste0(utr_exchange$receiving_tx[i], "longUTR"),
    transcript_source = "local",
    transcript_biotype = "extendedUTR",
    tag = "extendedUTR",
    transcript_support_level = 1,
    exon_number = seq_len(length(tmp)),
    exon_id = paste0(utr_exchange$receiving_tx[i], "longUTR_E", seq_len(length(tmp))),
    exon_version = 1, protein_id = NA_character_,
    protein_version = NA_character_,
    ccds_id = NA_character_
  )
  tmp
}), "GRangesList")

tx_expanded <- endoapply(ebt_expanded, function(x) {
  tmp <- range(x)
  mcols(tmp) <- data.frame(
    source = "local",
    type = "transcript",
    score = NA_real_,
    phase = NA_integer_,
    gene_id = x$gene_id[1],
    gene_version = x$gene_version[i],
    gene_name = x$gene_name[i],
    gene_source = x$gene_source[i],
    gene_biotype = x$gene_biotype[i],
    transcript_id = x$transcript_id[i],
    transcript_version = 1,
    transcript_name = x$transcript_id[i],
    transcript_source = "local",
    transcript_biotype = "extendedUTR",
    tag = "extendedUTR",
    transcript_support_level = 1,
    exon_number = NA_character_,
    exon_id = NA_character_,
    exon_version = 1, protein_id = NA_character_,
    protein_version = NA_character_,
    ccds_id = NA_character_
  )
  tmp
})

## ========================================================================== ##
## Add to original GRanges object
## ========================================================================== ##
## Exons
gtf_expanded <- c(gtf, unlist(ebt_expanded), unlist(tx_expanded))

## ========================================================================== ##
## Sort gtf
## ========================================================================== ##
gtf_expanded <- data.frame(gtf_expanded, stringsAsFactors = FALSE)
genelevels <- unique(gtf_expanded %>% 
  dplyr::arrange(as.numeric(seqnames), start) %>%
  dplyr::pull(gene_id))
txlevels <- c(NA_character_, setdiff(unique(gtf_expanded$transcript_id), NA))
typelevels <- c("gene", "transcript", "exon")
gtf_expanded <- gtf_expanded %>% 
  dplyr::mutate(gene_id = factor(gene_id, levels = genelevels),
                transcript_id = factor(transcript_id, levels = txlevels, exclude = NULL),
                type = factor(type, levels = typelevels)) %>%
  dplyr::arrange(gene_id, transcript_id, type)

gtfout <- GRanges(seqnames = gtf_expanded$seqnames,
                  ranges = IRanges(start = gtf_expanded$start,
                                   end = gtf_expanded$end),
                  strand = gtf_expanded$strand)
mcols(gtfout) <- gtf_expanded %>% dplyr::select(-seqnames, -start, -end, -width,
                                                -strand)
gtf_expanded <- gtfout

## ========================================================================== ##
## Save tx2gene and expanded gtf
## ========================================================================== ##
tx2gene <- as.data.frame(mcols(subset(gtf_expanded, 
                                      type == "transcript")))[, c("transcript_id", "gene_id")] %>%
  dplyr::distinct() %>%
  dplyr::rename(tx = transcript_id, gene = gene_id) %>%
  dplyr::mutate(symbol = gene)

rtracklayer::export(gtf_expanded, con = paste0(outbase, ".gtf"), format = "gtf")
saveRDS(tx2gene, file = paste0(outbase, "_tx2gene.rds"))

## ========================================================================== ##
## Get sequence of new transcripts
## ========================================================================== ##
genome <- BSgenome.Hsapiens.NCBI.GRCh38
ebt_expanded_seq <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, transcripts = ebt_expanded
)

## Add to original set of transcripts
ebt_expanded_seq <- ebt_expanded_seq[gsub("longUTR", "", names(ebt_expanded_seq)) %in% 
                                       names(txfasta)]
outfasta <- c(txfasta, ebt_expanded_seq)

writeXStringSet(outfasta, filepath = paste0(outbase, ".fasta"))

date()
sessionInfo()
