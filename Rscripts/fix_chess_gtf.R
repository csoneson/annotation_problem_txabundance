################################################################################
##                                                                            ##
## Fix CHESS gtf to be compatible with the analysis framework                 ##
##                                                                            ##
## Inputs:                                                                    ##
## * ingtf: input gtf from CHESS                                              ##
## * ingenes: list of genes present in CHESS gtf file, with Entrez IDs        ##
## * outgtf: output gtf file                                                  ##
## * outinfo: output file containing mapping between gene IDs and symbols     ##
##                                                                            ##
## Outputs:                                                                   ##
## * Modified gtf file                                                        ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(ingtf)
print(ingenes)
print(outgtf)
print(outinfo)

suppressPackageStartupMessages({
  library(GenomeInfoDb)
  library(rtracklayer)
  library(readr)
  library(org.Hs.eg.db)
  library(dplyr)
})

## Read list of genes and add Ensembl ID where it is available
genes <- read.delim(ingenes, header = TRUE, as.is = TRUE)
map <- AnnotationDbi::select(org.Hs.eg.db, keys = as.character(genes$RefSeq_GeneID), 
                             keytype = "ENTREZID", columns = c("ENTREZID", "ENSEMBL"))
genes$ensembl_id <- map$ENSEMBL[match(genes$RefSeq_GeneID, map$ENTREZID)]


## Read original gtf file
tmp <- read_tsv(ingtf, comment = "##", col_names = FALSE, col_types = "cccddcccc")
tmp <- as.data.frame(tmp)

## Convert all ranges to integers and write back to file
tmp$X4 <- as.integer(tmp$X4)
tmp$X5 <- as.integer(tmp$X5)
rn <- round(1e6 * runif(1))
write.table(tmp, file = paste0("tmp", rn, ".gff"), row.names = FALSE,
            col.names = FALSE, quote = FALSE, sep = "\t")

## Read back with rtracklayer and remove intermediate file
gtf <- import(paste0("tmp", rn, ".gff"), format = "gff")
unlink(paste0("tmp", rn, ".gff"))

## Convert sequence names to Ensembl format
seqlevelsStyle(gtf) <- "NCBI"

gtf$gene_id <- gtf$ID
gtf$gene_id[gtf$type %in% c("transcript", "tRNA", "rRNA")] <- 
  as.character(gtf$Parent[gtf$type %in% c("transcript", "tRNA", "rRNA")])
gtf$gene_id[gtf$type %in% c("exon", "CDS")] <- 
  as.character(gtf$Parent[as.numeric(match(gtf$Parent[gtf$type %in% c("exon", "CDS")], gtf$ID))])

gtf$gene_name <- genes$ensembl_id[match(gtf$gene_id, genes$GFF_ID)]

gtf$transcript_name <- gtf$transcript_id

gtf$transcript_id <- NA_character_
gtf$transcript_id[gtf$type %in% c("transcript", "tRNA", "rRNA")] <- 
  as.character(gtf$ID[gtf$type %in% c("transcript", "tRNA", "rRNA")])
gtf$transcript_id[gtf$type %in% c("exon", "CDS")] <- 
  as.character(gtf$Parent[gtf$type %in% c("exon", "CDS")])

gtf$exon_id <- NA_character_
gtf$exon_id[gtf$type == "exon"] <- paste0("E", seq_len(length(which(gtf$type == "exon"))))

mcols(gtf) <- mcols(gtf)[, c("source", "type", "score", "phase", "gene_id", 
                             "transcript_id", "exon_id", "gene_name", "transcript_name")]

gtf <- subset(gtf, source != "tRNAscan-SE")

## Save a file with the conversion between gene ID and "symbol" (here, Ensembl ID)
saveRDS(as.data.frame(mcols(gtf)) %>% dplyr::select(gene_id, gene_name) %>%
          dplyr::rename(gene = gene_id, symbol = gene_name),
        file = outinfo)

rtracklayer::export(gtf, outgtf, format = "gtf")

date()
sessionInfo()
