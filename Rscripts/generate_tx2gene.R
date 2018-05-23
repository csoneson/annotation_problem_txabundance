################################################################################
##                                                                            ##
## Generate transcript-to-gene conversion table                               ##
##                                                                            ##
## Inputs:                                                                    ##
## * cdna: cDNA fasta file (from Ensembl)                                     ##
## * ncrna: ncRNA fasta file (from Ensembl)                                   ##
## * cds: CDS fasta file (from Ensembl)                                       ##
## * outrds: output file                                                      ##
##                                                                            ##
## Outputs:                                                                   ##
## * A transcript-to-gene conversion table                                    ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(cdna)
print(ncrna)
print(cds)
print(outrds)

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
})

cdna <- readDNAStringSet(cdna)
ncrna <- readDNAStringSet(ncrna)
cds <- readDNAStringSet(cds)

tx <- c(cdna, ncrna, cds)

tx2gene <- data.frame(t(sapply(as.character(names(tx)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  symbol <- gsub("^gene_symbol:", "", a[grep("^gene_symbol:", a)])
  gene_biotype <- gsub("^gene_biotype:", "", a[grep("^gene_biotype:", a)])
  tx_biotype <- gsub("^transcript_biotype:", "", a[grep("^transcript_biotype:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), 
    gene = ifelse(length(gene) != 0, gene, NA),
    symbol = ifelse(length(symbol) != 0, symbol, NA),
    gene_biotype = ifelse(length(gene_biotype) != 0, gene_biotype, NA),
    tx_biotype = ifelse(length(tx_biotype) != 0, tx_biotype, NA))
})), stringsAsFactors = FALSE)
rownames(tx2gene) <- NULL
tx2gene <- distinct(tx2gene)

saveRDS(tx2gene, file = outrds)

sessionInfo()
date()
