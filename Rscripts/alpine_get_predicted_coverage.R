args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(bam)  ## bam file
print(biasmodels)  ## bias model object (output of alpine_fitbiasmodel.R)
print(ncores)  ## number of cores for parallel computations
print(outrds)  ## output file

suppressPackageStartupMessages(library(alpine))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))
suppressPackageStartupMessages(library(parallel))

## Read bias model parameters and exon-by-transcript objects
biasmodels <- readRDS(biasmodels)
fitpar <- biasmodels$fitpar
ebt0 <- biasmodels$ebt0

## Estimate average fragment length
avefraglength <- sum(fitpar$`1`$fraglen.density$x * fitpar$`1`$fraglen.density$y/
                       sum(fitpar$`1`$fraglen.density$y))

## Load bam file 
bam.files <- bam
names(bam.files) <- "1"

## Get all transcripts
transcripts <- names(ebt0)
names(transcripts) <- transcripts

## Predict coverage for each transcript
res <- mclapply(transcripts, function(tx) {
  message(tx)
  ## Get transcript model
  txmod <- ebt0[[tx]]
  txmods <- sort(ebt0[[tx]])
  
  m <- tryCatch({
    a <- predictCoverage(gene = txmod,
                         bam.files = bam.files,
                         fitpar = fitpar,
                         genome = Hsapiens,
                         model.names = "all")
    a$`1`$pred.cov$all
  }, 
  error = function(e) {
    ## Assume uniform coverage
    Rle(rep(1, sum(width(txmod)) - 1))
  })
  
  junctions <- GenomicRanges::setdiff(range(txmods), txmods)
  if (all(strand(txmods) == "+")) {
    junctionpos <- cumsum(width(txmods))
    junctionpos <- junctionpos[-length(junctionpos)]
    mcols(junctions)$pred.cov <- as.numeric(m)[junctionpos]
    strand <- "+"
  } else if (all(strand(txmods) == "-")) {
    junctionpos <- cumsum(width(rev(txmods)))
    junctionpos <- junctionpos[-length(junctionpos)]
    mcols(junctions)$pred.cov <- rev(as.numeric(m)[junctionpos])
    strand <- "-"
  } else {
    strand <- "mixed"
  }
  
  list(pred.cov = m, strand = strand, junctions = junctions, 
       avefraglength = avefraglength)
  
}, mc.preschedule = FALSE, mc.cores = ncores)

saveRDS(res, file = outrds)

sessionInfo()
date()
