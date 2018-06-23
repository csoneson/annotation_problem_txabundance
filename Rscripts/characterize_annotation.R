################################################################################
##                                                                            ##
## Characterize annotation                                                    ##
##                                                                            ##
## Inputs:                                                                    ##
## * gtf: gtf file                                                            ##
## * txome: transcriptome fasta file                                          ##
## * outtxt: output file                                                      ##
##                                                                            ##
## Outputs:                                                                   ##
## * A text file with annotation characteristics                              ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(txome)
print(outtxt)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(dplyr)
  library(GenomicFeatures)
  library(GenomicRanges)
})

gtf <- rtracklayer::import(gtf, format = "gtf")
txome <- readDNAStringSet(txome)
tx2gene <- as.data.frame(gtf) %>% dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::distinct()

## Number of genes in gtf file
nbrGenesGtf <- length(unique(gtf$gene_id))

## Number of transcripts in gtf file
nbrTxGtf <- length(unique(gtf$transcript_id))

## Number of transcripts/gene in gtf file
gtftx <- subset(gtf, type == "transcript")
nbrTxPerGeneGtf <- table(gtftx$gene_id)
minNbrTxPerGeneGtf <- min(nbrTxPerGeneGtf)
meanNbrTxPerGeneGtf <- mean(nbrTxPerGeneGtf)
maxNbrTxPerGeneGtf <- max(nbrTxPerGeneGtf)
medianNbrTxPerGeneGtf <- median(nbrTxPerGeneGtf)
nbrNbrTxPerGeneGtf1 <- length(which(nbrTxPerGeneGtf == 1))

## Number of exons/transcript in gtf file
gtfex <- subset(gtf, type == "exon")
nbrExonPerTxGtf <- table(gtfex$transcript_id)
minNbrExonPerTxGtf <- min(nbrExonPerTxGtf)
meanNbrExonPerTxGtf <- mean(nbrExonPerTxGtf)
maxNbrExonPerTxGtf <- max(nbrExonPerTxGtf)
medianNbrExonPerTxGtf <- median(nbrExonPerTxGtf)
nbrNbrExonPerTxGtf1 <- length(which(nbrExonPerTxGtf == 1))

## Transcript length
txLengthsGtf <- as.data.frame(gtfex) %>% dplyr::group_by(transcript_id) %>%
  dplyr::summarize(length = sum(width)) %>% dplyr::pull(length)
minTxLengthGtf <- min(txLengthsGtf)
meanTxLengthGtf <- mean(txLengthsGtf)
maxTxLengthGtf <- max(txLengthsGtf)
medianTxLengthGtf <- median(txLengthsGtf)

## Number of junctions per gene
txdb <- makeTxDbFromGRanges(gtf)
ebt <- exonsBy(txdb, by = "tx", use.names = TRUE)
jbt <- endoapply(ebt, function(w) GenomicRanges::setdiff(range(w), w))
junctions <- as.data.frame(jbt) %>%
  dplyr::left_join(tx2gene, by = c("group_name" = "transcript_id"))
junctionsByGeneGtf <- junctions %>% dplyr::select(gene_id, seqnames, start, end, strand) %>%
  dplyr::distinct()
nbrJunctionsPerGeneGtf <- table(junctionsByGeneGtf$gene_id)
minNbrJunctionsPerGeneGtf <- min(nbrJunctionsPerGeneGtf)
meanNbrJunctionsPerGeneGtf <- mean(nbrJunctionsPerGeneGtf)
maxNbrJunctionsPerGeneGtf <- max(nbrJunctionsPerGeneGtf)
medianNbrJunctionsPerGeneGtf <- median(nbrJunctionsPerGeneGtf)
nbrNbrJunctionsPerGeneGtf1 <- sum(nbrJunctionsPerGeneGtf == 1)
nbrNbrJunctionsPerGeneGtf0 <- length(setdiff(unique(gtf$gene_id), names(nbrJunctionsPerGeneGtf)))

## Number of transcripts not present in the gtf
txsfa <- sapply(strsplit(names(txome), " "), .subset, 1)
idx <- grep("^STRG\\.|^CHS\\.", txsfa, invert = TRUE)
txsfa[idx] <- gsub("\\.[0-9]+$", "", txsfa[idx])
nbrTxNotInGtf <- length(setdiff(txsfa, gtf$transcript_id))

write.table(t(data.frame(nbrGenesGtf = nbrGenesGtf,
                         nbrTxGtf = nbrTxGtf,
                         minNbrTxPerGeneGtf = minNbrTxPerGeneGtf,
                         meanNbrTxPerGeneGtf = meanNbrTxPerGeneGtf,
                         maxNbrTxPerGeneGtf = maxNbrTxPerGeneGtf,
                         medianNbrTxPerGeneGtf = medianNbrTxPerGeneGtf,
                         nbrNbrTxPerGeneGtf1 = nbrNbrTxPerGeneGtf1,
                         minNbrExonPerTxGtf = minNbrExonPerTxGtf,
                         meanNbrExonPerTxGtf = meanNbrExonPerTxGtf,
                         maxNbrExonPerTxGtf = maxNbrExonPerTxGtf,
                         medianNbrExonPerTxGtf = medianNbrExonPerTxGtf,
                         nbrNbrExonPerTxGtf1 = nbrNbrExonPerTxGtf1,
                         minTxLengthGtf = minTxLengthGtf,
                         meanTxLengthGtf = meanTxLengthGtf,
                         maxTxLengthGtf = maxTxLengthGtf,
                         medianTxLengthGtf = medianTxLengthGtf,
                         minNbrJunctionsPerGeneGtf = minNbrJunctionsPerGeneGtf,
                         meanNbrJunctionsPerGeneGtf = meanNbrJunctionsPerGeneGtf,
                         maxNbrJunctionsPerGeneGtf = maxNbrJunctionsPerGeneGtf,
                         medianNbrJunctionsPerGeneGtf = medianNbrJunctionsPerGeneGtf,
                         nbrNbrJunctionsPerGeneGtf1 = nbrNbrJunctionsPerGeneGtf1,
                         nbrNbrJunctionsPerGeneGtf0 = nbrNbrJunctionsPerGeneGtf0,
                         nbrTxNotInGtf = nbrTxNotInGtf,
                         stringsAsFactors = FALSE)),
            file = outtxt, row.names = TRUE, col.names = FALSE,
            quote = FALSE, sep = "\t")

date()
sessionInfo()
