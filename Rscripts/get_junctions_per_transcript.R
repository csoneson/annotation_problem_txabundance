args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(outrds)

suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicRanges))

## Construct ensembldb object from gtf file
ensDbFromGtf(gtf, path = dirname(outrds))

## Get list of transcripts
txdb <- EnsDb(paste0(dirname(outrds), "/", gsub("gtf", "sqlite", basename(gtf))))

## Get list of exons by transcript
ebt <- exonsBy(txdb, by = "tx")

## Junctions by transcript
jbt <- endoapply(ebt, function(w) {
  setdiff(range(w), w)
})

saveRDS(jbt, file = outrds)

sessionInfo()
date()