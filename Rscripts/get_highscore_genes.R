args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(scorerds)
print(outtxt)

suppressPackageStartupMessages({
  library(dplyr)
})

scores <- readRDS(scorerds)$genes

scores_sub <- 
  scores %>% dplyr::filter(count > 1000 & intron_exon_ratio < 0.1 & 
                             uniqjuncfraction > 0.9 & uniqjuncreads > 25) %>%
  dplyr::group_by(method) %>% dplyr::arrange(desc(score)) %>% dplyr::top_n(50) %>%
  dplyr::ungroup() %>% dplyr::group_by(gene) %>% dplyr::mutate(meanscore = mean(score)) %>% 
  dplyr::arrange(desc(meanscore), method) %>%
  dplyr::select(gene, method, count, nbr_transcripts, 
                length_diff_3putrs_samestart, intron_exon_ratio,
                uniqjuncreads, mmjuncreads, 
                uniqjuncfraction, nbr_junctions_in_gene, score)

write.table(scores_sub, file = outtxt, row.names = FALSE, col.names = TRUE,
            sep = "\t", quote = FALSE)

## Shared high-scoring genes
scores_shared <- 
  scores %>% dplyr::filter(method %in% c("hera", "kallisto", "RSEM", "Salmon",
                                         "SalmonBWA", "StringTie") & 
                             count > 1000 & score >= 0.6 & uniqjuncreads > 25 & 
                             uniqjuncfraction > 0.75 & intron_exon_ratio < 0.1) %>% 
  dplyr::group_by(gene) %>% dplyr::mutate(n = length(method)) %>% 
  dplyr::ungroup() %>% dplyr::filter(n == length(unique(method))) %>%
  dplyr::select(gene, method, count, nbr_transcripts, 
                length_diff_3putrs_samestart, intron_exon_ratio,
                uniqjuncreads, mmjuncreads, 
                uniqjuncfraction, nbr_junctions_in_gene, score)

print(length(unique(scores_shared$gene)))

write.table(scores_shared, file = gsub("\\.txt", "_shared.txt", outtxt), 
            row.names = FALSE, col.names = TRUE,
            sep = "\t", quote = FALSE)

date()
sessionInfo()
