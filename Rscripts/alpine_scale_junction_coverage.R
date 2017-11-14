args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(predcovrds)  ## predicted transcript and junction coverage
print(txquants)  ## transcript quantification file
print(quantreadscript)  ## script containing a function "read_quant()"
print(strandspec)  ## "yes" or "no"
print(outrds)  ## output file

suppressPackageStartupMessages(library(alpine))
suppressPackageStartupMessages(library(parallel))

source(quantreadscript)

predcovs <- readRDS(predcovrds)
quants <- read_quant(file = txquants, avefraglength = predcovs[[1]]$avefraglength)

transcripts <- names(predcovs)
names(transcripts) <- transcripts

scaledcovs <- lapply(transcripts, function(tx) {
  ab <- quants$count[quants$transcript == tx]
  m <- predcov[[tx]]s$junctions
  m$pred.cov <- m$pred.cov / sum(predcovs[[tx]]$pred.cov) * 
    ab * predcovs[[tx]]$avefraglength
  as.data.frame(m)
})

allcovs <- do.call(rbind, scaledcovs)
if (strandspec == "no") {
  allcovs$strand <- "*"
}
allcovs <- allcovs %>%
  dplyr::group_by(seqnames, start, end, strand) %>%
  dplyr::summarize(pred.cov = sum(pred.cov))

saveRDS(list(scaledcovs = scaledcovs, allcovs = allcovs, quants = quants), 
        file = outrds)

sessionInfo()
date()