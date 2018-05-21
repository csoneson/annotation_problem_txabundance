################################################################################
##                                                                            ##
## Fit bias model with alpine                                                 ##
##                                                                            ##
## Inputs:                                                                    ##
## * gtf: gtf file                                                            ##
## * bam: bam file, used to extract "medium to highly expressed genes" to fit ##
##        the bias model                                                      ##
## * readlength: read length                                                  ##
## * minsize, maxsize: minimal and maximal fragment length                    ##
## * organism: organism                                                       ##
## * genomeVersion: genome version                                            ##
## * version: Ensembl version                                                 ##
## * outdir: output directory                                                 ##
## * subsample: whether to use only a subset of the "eligible" genes to fit   ##
##              the model                                                     ##
##                                                                            ##
## Outputs:                                                                   ##
## * A bias model and a list of exons by transcript                           ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(bam)  ## used to extract "medium to highly expressed genes" to fit the bias model
print(readlength)
print(minsize)
print(maxsize)
print(organism)
print(genomeVersion)
print(version)
print(outdir)
print(subsample)

suppressPackageStartupMessages({
  library(ensembldb)
  library(alpine)
  library(Rsamtools)
  library(GenomicFeatures)
  library(GenomicAlignments)
  library(BSgenome.Hsapiens.NCBI.GRCh38)
  library(dplyr)
})

## Construct ensembldb object from gtf file and get TxDb
txdb <- 
  tryCatch({
    ensDbFromGtf(gtf, path = outdir, organism = organism, 
                 genomeVersion = genomeVersion, version = version)
    EnsDb(paste0(outdir, "/", gsub("gtf", "sqlite", basename(gtf))))
  }, error = function(e) {
    GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf", organism = gsub("_", " ", organism))
  })
txdb

## Get list of transcripts
if (class(txdb) == "EnsDb") {
  txdf <- transcripts(txdb, return.type = "DataFrame")  ## data frame format
  txps <- transcripts(txdb)  ## GRanges format
} else {
  txps <- GenomicFeatures::transcripts(txdb, columns = c("tx_name", "gene_id"), use.names = TRUE)
  txps$gene_id <- unlist(txps$gene_id)
  txdf <- DataFrame(tx_id = as.character(txps$tx_name), gene_id = as.character(unlist(txps$gene_id)),
                    tx_seq_start = as.integer(start(txps)), tx_seq_end = as.integer(end(txps)),
                    tx_name = as.character(txps$tx_name))
}
txdf
txps

## Select genes with a single isoform
tab <- table(txdf$gene_id)
one.iso.genes <- names(tab)[tab == 1]
length(one.iso.genes)
length(tab)

## Get list of exons by transcript
if (class(txdb) == "EnsDb") {
  (ebt0 <- exonsBy(txdb, by = "tx"))
} else { 
  (ebt0 <- exonsBy(txdb, by = "tx", use.names = TRUE))
}

## Get transcript names for genes with a single isoform
one.iso.txs <- txdf$tx_id[txdf$gene_id %in% one.iso.genes]

## Extract these transcripts for use in fitting bias model
ebt.fit <- ebt0[one.iso.txs]
ebt.fit <- keepStandardChromosomes(ebt.fit, pruning.mode = "coarse")

## Filter small genes and long genes
min.bp <- 600
max.bp <- 7000
gene.lengths <- sum(width(ebt.fit))
summary(gene.lengths)

ebt.fit <- ebt.fit[gene.lengths > min.bp & gene.lengths < max.bp]
length(ebt.fit)

## Sample 1000 genes to use for the fitting of the bias model
if (subsample) {
  set.seed(1)
  ebt.fit <- ebt.fit[sample(length(ebt.fit), 1000)]
}

## Read bam file
bam.files <- bam
names(bam.files) <- "1"

## Subset to medium-to-high expressed genes
txps.fit <- sort(txps[names(ebt.fit)])
cts <- countBam(BamFile(bam.files), param = ScanBamParam(which = txps.fit))
mcols(txps.fit)$cts <- cts$records
txps.fit <- txps.fit[names(ebt.fit)]
sum(txps.fit$cts > 500 & txps.fit$cts < 10000)
ebt.fit <- ebt.fit[txps.fit$cts > 500 & txps.fit$cts < 10000]
length(ebt.fit)

## Get fragment width and read length
m <- sort(which(sapply(ebt.fit, length) > 1))
m0 <- 0
w <- NULL
while(is.null(w)) {
  w <- tryCatch({
    m0 <- m0 + 1
    getFragmentWidths(bam.files, ebt.fit[[m[m0]]])
  }, error = function(e) {
    NULL
  })
}
quantile(w, c(.025, .975))
getReadLength(bam.files)

## Names of genes to retain
gene.names <- names(ebt.fit)
names(gene.names) <- gene.names

## Build fragment types for selected genes
system.time({
  fragtypes <- lapply(gene.names, function(gene.name) {
    buildFragtypes(exons = ebt.fit[[gene.name]],
                   genome = Hsapiens,
                   readlength = readlength,
                   minsize = minsize,
                   maxsize = maxsize,
                   gc.str = FALSE)
  })
})
object.size(fragtypes)/1e6

## Define models to fit
models <- list(
  "all" = list(
    formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + 
    ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) + gene",
    offset = c("fraglen", "vlmm")
  )
)

## Fit bias model
system.time({
  fitpar <- lapply(bam.files, function(bf) {
    fitBiasModels(genes = ebt.fit,
                  bam.file = bf,
                  fragtypes = fragtypes,
                  genome = Hsapiens,
                  models = models,
                  readlength = readlength,
                  minsize = minsize,
                  maxsize = maxsize)
  })
})

## Diagnostic plots
pdf(paste0(outdir, "/alpine_fitbiasmodel_plots.pdf"))
plotFragLen(fitpar)
plotGC(fitpar, model = "all")
dev.off()

## Save bias model
saveRDS(list(fitpar = fitpar, ebt0 = ebt0, txps = txps), 
        file = paste0(outdir, "/alpine_fitbiasmodel.rds"))

date()
sessionInfo()
