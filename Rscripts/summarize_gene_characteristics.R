args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(quantsf)
print(gtf)
print(tx2gene)
print(outrds)

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(dplyr))

tx2gene <- readRDS(tx2gene)

## Read transcript abundances
quantsf <- read.delim(quantsf, header = TRUE, as.is = TRUE)

## Add gene ID
quantsf$gene <- tx2gene$gene[match(quantsf$Name, tx2gene$tx)]

## Add number of junctions per transcript
gtf <- import(gtf)
gtf <- subset(gtf, type == "exon")
nbrexon <- as.data.frame(gtf) %>% group_by(transcript_id) %>% 
  dplyr::summarize(nex = length(exon_id),
                   totlength = sum(width))

## Summarize by gene
genechar <- quantsf %>% dplyr::mutate(Name = gsub("\\.[0-9]+", "", Name),
                                      gene = gsub("\\.[0-9]+", "", gene)) %>%
  dplyr::left_join(nbrexon, by = c("Name" = "transcript_id")) %>%
  group_by(gene) %>% dplyr::mutate(TPMrel = TPM/sum(TPM)) %>% 
  dplyr::mutate(TPMrel = replace(TPMrel, is.na(TPMrel), 0)) %>% 
  dplyr::summarize(ntx = length(Name),
                   TPM = sum(TPM),
                   count = sum(NumReads),
                   ntxabove5p = sum(TPMrel > 0.05),
                   maxnbrex = max(nex),
                   maxlength = max(totlength))

saveRDS(genechar %>% dplyr::filter(maxnbrex > 1 & 
                                     count > 1000 & 
                                     maxlength > 350),
        file = outrds)

sessionInfo()
date()

