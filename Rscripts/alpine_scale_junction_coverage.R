################################################################################
##                                                                            ##
## Scale predicted junction coverages based on estimated transcript           ##
## abundances                                                                 ##
##                                                                            ##
## Inputs:                                                                    ##
## * predcovrds: predicted transcript coverage patterns and junction          ##
##               coverages (output from alpine_get_predicted_coverage.R)      ##
## * txquants: file with transcript abundances                                ##
## * quantreadscript: script to read the transcript quantifications. Should   ##
##                    contain a function named read_quant()                   ##
## * strandspec: yes or no, is the data strand-specific                       ##
## * tx2gene: data frame with transcript-to-gene conversion info              ##
## * method: method name to add to the quantification table                   ##
## * outrds: output file                                                      ##
##                                                                            ##
## Outputs:                                                                   ##
## * Predicted junction coverages for all junctions, and a table with all     ##
##   transcript abundances                                                    ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(predcovrds)
print(txquants)
print(quantreadscript)
print(strandspec)
print(tx2gene)
print(method)
print(outrds)

suppressPackageStartupMessages({
  library(alpine)
  library(parallel)
  library(dplyr)
})

source(quantreadscript)

predcovs <- readRDS(predcovrds)
quants <- read_quant(file = txquants, avefraglength = predcovs[[1]]$avefraglength)
tx2gene <- readRDS(tx2gene)

idxtx <- grep("^STRG\\.", tx2gene$tx, invert = TRUE)
tx2gene$tx[idxtx] <- gsub("\\.[0-9]+$", "", tx2gene$tx[idxtx])

idxgene <- grep("^STRG\\.", tx2gene$gene, invert = TRUE)
tx2gene$gene[idxgene] <- gsub("\\.[0-9]+$", "", tx2gene$gene[idxgene])

transcripts <- names(predcovs)
names(transcripts) <- transcripts

## Go through all transcripts and scale predicted junction coverage by the
## estimated abundance
scaledcovs <- lapply(transcripts, function(tx) {
  tryCatch({
    ab <- quants$count[quants$transcript == tx]
    if (length(ab) == 0 || is.na(ab)) ab <- 0  ## if the transcript is not present in the quantification file
    m <- predcovs[[tx]]$junctions
    m$pred.cov <- m$pred.cov / max(1e-10, sum(predcovs[[tx]]$pred.cov)) * 
      ab * predcovs[[tx]]$avefraglength
    as.data.frame(m) %>% dplyr::mutate(transcript = tx)
  }, error = function(e) NULL)
})

## Combine estimates for all transcripts
allcovs <- do.call(rbind, scaledcovs)
if (strandspec == "no") {
  allcovs$strand <- "*"
}

## Sum coverages of the same junction from different transcripts. Each junction
## is present once per gene it is included in.
allcovs <- allcovs %>%
  dplyr::group_by(seqnames, start, end, strand) %>%
  dplyr::mutate(pred.cov = sum(pred.cov)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(tx2gene %>% dplyr::select(tx, gene, symbol), by = c("transcript" = "tx")) %>%
  dplyr::group_by(seqnames, start, end, strand, gene) %>%
  dplyr::mutate(transcript = paste(transcript, collapse = ",")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>% 
  dplyr::mutate(method = method)

## Get coverage "note" for each transcript
covnotes <- sapply(transcripts, function(tx) {
  predcovs[[tx]]$note
})

## Add gene info to transcript quant table
quants <- quants %>% dplyr::left_join(tx2gene %>% dplyr::select(tx, gene, symbol), 
                                      by = c("transcript" = "tx")) %>%
  dplyr::mutate(method = method) %>% 
  dplyr::left_join(data.frame(transcript = names(covnotes), 
                              covnote = covnotes, 
                              stringsAsFactors = FALSE), by = "transcript") %>%
  dplyr::select(transcript, gene, symbol, count, TPM, method, covnote)

saveRDS(list(scaledcovs = scaledcovs, allcovs = allcovs, quants = quants), 
        file = outrds)

sessionInfo()
date()