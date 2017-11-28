args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(combcovrds)
print(outrds)

suppressPackageStartupMessages(library(dplyr))

x <- readRDS(combcovrds)

quants <- x$allquants

## Summarize on gene level
quants_gene <- quants %>% dplyr::group_by(gene, method) %>%
  dplyr::mutate(TPMrel = TPM/sum(TPM)) %>%
  dplyr::mutate(TPMrel = replace(TPMrel, is.na(TPMrel), 0)) %>%
  dplyr::summarize(count = sum(count),
                   TPM = sum(TPM),
                   nbr_expressed_transcripts = sum(TPM > 0),
                   nbr_expressed_transcripts_5p = sum(TPMrel > 0.05)) %>%
  dplyr::ungroup()

saveRDS(quants_gene, file = outrds)

sessionInfo()
date()
