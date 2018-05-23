################################################################################
##                                                                            ##
## Filter a StringTie output gtf to exclude transcripts with unknown strand   ##
##                                                                            ##
## Inputs:                                                                    ##
## * ingtf: input gtf (from StringTie)                                        ##
## * outgtf: output gtf                                                       ##
##                                                                            ##
## Outputs:                                                                   ##
## * A gtf file where transcripts with unknown strand have been excluded, and ##
##   gene_id and transcript_ids have been replaced with ref_gene_id and       ##
##   reference_id whenever applicable.                                        ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(ingtf)
print(outgtf)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(readr)
})

## Read the file with read_tsv and save as a temporary file without the top
## commented lines, that can be read with rtracklayer::import
tmp <- readr::read_tsv(ingtf, comment = "#", col_names = FALSE, col_types = "ccccccccc")
rn <- round(1e6*runif(1))
write.table(tmp, file = paste0("tmp", rn, ".gtf"), col.names = FALSE, row.names = FALSE, 
            sep = "\t", quote = FALSE)
gtf <- import(paste0("tmp", rn, ".gtf"), format = "gtf")
unlink(paste0("tmp", rn, ".gtf"))

## Filter out transcripts with unknown strand
gtf <- gtf[strand(gtf) != "*"]

## Replace gene_id and transcript_id with ref_gene_id and reference_id when possible
nona <- !is.na(gtf$reference_id)
gtf$transcript_id[nona] <- gtf$reference_id[nona]
nona <- !is.na(gtf$ref_gene_id)
gtf$gene_id[nona] <- gtf$ref_gene_id[nona]

## Write filtered gtf
export(gtf, outgtf, format = "gtf")

date()
sessionInfo()
