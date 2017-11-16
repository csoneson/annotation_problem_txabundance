args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(predcovrds)  ## predicted transcript and junction coverage
print(txquants)  ## transcript quantification file
print(quantreadscript)  ## script containing a function "read_quant()"
print(strandspec)  ## "yes" or "no"
print(tx2gene)  ## link transcripts to genes
print(method)  ## method ID to add to the quantification table
print(outrds)  ## output file

suppressPackageStartupMessages(library(alpine))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))

source(quantreadscript)

predcovs <- readRDS(predcovrds)
quants <- read_quant(file = txquants, avefraglength = predcovs[[1]]$avefraglength)
tx2gene <- readRDS(tx2gene)
tx2gene$tx <- gsub("\\.[0-9]+", "", tx2gene$tx)
tx2gene$gene <- gsub("\\.[0-9]+", "", tx2gene$gene)

transcripts <- names(predcovs)
names(transcripts) <- transcripts

scaledcovs <- lapply(transcripts, function(tx) {
  tryCatch({
    ab <- quants$count[quants$transcript == tx]
    m <- predcovs[[tx]]$junctions
    m$pred.cov <- m$pred.cov / sum(predcovs[[tx]]$pred.cov) * 
      ab * predcovs[[tx]]$avefraglength
    as.data.frame(m) %>% dplyr::mutate(transcript = tx)
  }, error = function(e) NULL)
})

allcovs <- do.call(rbind, scaledcovs)
if (strandspec == "no") {
  allcovs$strand <- "*"
}
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

saveRDS(list(scaledcovs = scaledcovs, allcovs = allcovs, quants = quants), 
        file = outrds)

sessionInfo()
date()