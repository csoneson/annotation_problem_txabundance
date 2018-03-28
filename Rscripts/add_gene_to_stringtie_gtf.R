args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Add the "gene" rows to a StringTie gtf

print(ingtf)
print(outgtf)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
})

gtf <- import(ingtf)
gtf <- as.data.frame(gtf)

gtfgene <- gtf %>% dplyr::group_by(gene_id) %>%
  dplyr::summarize(seqnames = seqnames[1], 
                   start = min(start),
                   end = max(end),
                   width = max(end) - min(start) + 1,
                   strand = strand[1],
                   source = "StringTie",
                   type = "gene",
                   score = 1000,
                   phase = NA,
                   transcript_id = NA_character_,
                   cov = NA_character_,
                   FPKM = NA_character_,
                   TPM = NA_character_,
                   exon_number = NA_character_,
                   reference_id = NA_character_,
                   ref_gene_id = NA_character_,
                   ref_gene_name = NA_character_) %>%
  as.data.frame() %>%
  dplyr::select(colnames(gtf))

gtfcombined <- rbind(gtf, gtfgene)

gtfout <- GRanges(seqnames = gtfcombined$seqnames,
                  ranges = IRanges(start = gtfcombined$start,
                                   end = gtfcombined$end),
                  strand = gtfcombined$strand)
mcols(gtfout) <- gtfcombined %>% dplyr::select(-seqnames, -start, -end, -width,
                                               -strand)
gtfout <- sort(gtfout)

rtracklayer::export(gtfout, outgtf, format = "gtf")

sessionInfo()
date()
