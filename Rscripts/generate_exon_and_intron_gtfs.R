args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(ingtf)
print(outexongtf)
print(outintrongtf)

suppressPackageStartupMessages({
  library(rtracklayer)
})

ingtf <- rtracklayer::import(ingtf)

## Keep only exons
ingtf <- subset(ingtf, type == "exon")

## Split by gene
ingtf <- split(ingtf, ingtf$gene_id)

## For each gene, get the full range as well as the flattened exons
rangegtf <- endoapply(ingtf, function(w) {
  g <- w$gene_id[1]
  w <- range(w)
  mcols(w) <- data.frame(gene_id = g, stringsAsFactors = FALSE)
  w
})

## Reduce the exons for each gene
flatexongtf <- unlist(endoapply(ingtf, function(w) {
  g <- w$gene_id[1]
  w <- reduce(w)
  mcols(w) <- data.frame(gene_id = g, exon_id = paste0("e", seq_along(w)),
                         type = "exon",
                         stringsAsFactors = FALSE)
  w
}))

## For each gene, get introns by subtracting all exons from the gene range
introngtf <- unlist(endoapply(rangegtf, function(w) {
  g <- w$gene_id[1]
  w <- GenomicRanges::setdiff(w, flatexongtf)
  mcols(w) <- data.frame(gene_id = g, exon_id = paste0("e", seq_along(w)),
                         type = "exon",
                         stringsAsFactors = FALSE)
  w
}))

rtracklayer::export(flatexongtf, outexongtf, format = "gtf")
rtracklayer::export(introngtf, outintrongtf, format = "gtf")

date()
sessionInfo()
