args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(tx2gene)
print(gtf)
print(outrds)

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(dplyr))

## Extend tx2gene with the transcripts only included in the gtf file

tx2gene <- readRDS(tx2gene)
gtf <- import(gtf)

gtfsub <- as.data.frame(gtf) %>% dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, transcript_version, gene_id, gene_version, 
               gene_name, gene_biotype, transcript_biotype) %>%
  dplyr::mutate(transcript_id = paste0(transcript_id, ".", transcript_version),
                gene_id = paste0(gene_id, ".", gene_version)) %>%
  dplyr::select(-transcript_version, -gene_version) %>%
  dplyr::rename(tx_biotype = transcript_biotype,
                tx = transcript_id, 
                gene = gene_id,
                symbol = gene_name) %>%
  dplyr::select(tx, gene, symbol, gene_biotype, tx_biotype)

tx2gene <- tx2gene %>% dplyr::select(tx, gene, symbol, gene_biotype, tx_biotype)

tx2gene <- rbind(tx2gene, gtfsub) %>% dplyr::distinct()

saveRDS(tx2gene, file = outrds)

sessionInfo()
date()
