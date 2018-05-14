args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(scorerds)
print(targetmethod) ## Method to which all others will be compared
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

scores <- readRDS(scorerds)
scores <- scores$genes

## Consider only genes with uniqjuncreads > 25
tmpscore <- scores %>% dplyr::filter(uniqjuncreads > 25) %>% 
  dplyr::select(gene, method, score)
scoresub <- tmpscore %>%
  tidyr::spread(method, score) %>% as.data.frame() %>% 
  tibble::column_to_rownames("gene")

png(gsub("rds$", "png", outrds), width = 8, height = 8, unit = "in", res = 300)
plots <- lapply(setdiff(colnames(scoresub), targetmethod), function(cn) {
  ggplot(scoresub, aes_string(x = cn, y = targetmethod)) + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_point(alpha = 0.3, size = 0.3) + theme_bw()
})
print(cowplot::plot_grid(plotlist = plots))
dev.off()

scorecomp <- tmpscore %>%
  dplyr::left_join(tmpscore %>% dplyr::filter(method == targetmethod) %>%
                     dplyr::select(gene, score) %>%
                     dplyr::rename(targetscore = score)) %>%
  dplyr::filter(method != targetmethod) %>%
  dplyr::group_by(method) %>%
  dplyr::summarize(ImNmorethanp5 = sum(score - targetscore >= 0.5, na.rm = TRUE),
                   ImNbtw0andp5 = sum(score - targetscore > 0 & score - targetscore < 0.5, 
                                     na.rm = TRUE),
                   ImN0 = sum(score == targetscore, na.rm = TRUE),
                   NmIbtw0andp5 = sum(score - targetscore < 0 & score - targetscore > (-0.5), 
                                     na.rm = TRUE),
                   NmImorethanp5 = sum(score - targetscore <= (-0.5), na.rm = TRUE)) %>%
  tidyr::gather(category, nbr_genes, -method) %>%
  dplyr::mutate(category = replace(category, category == "ImN0", "score = Nanopore score"),
                category = replace(category, category == "ImNbtw0andp5", "0 < score - Nanopore score < 0.5"),
                category = replace(category, category == "ImNmorethanp5", "score - Nanopore score >= 0.5"),
                category = replace(category, category == "NmIbtw0andp5", "0 < Nanopore score - score < 0.5"),
                category = replace(category, category == "NmImorethanp5", "Nanopore score - score >= 0.5")) %>%
  dplyr::mutate(category = factor(category, levels = c("Nanopore score - score >= 0.5", 
                                                       "0 < Nanopore score - score < 0.5",
                                                       "score = Nanopore score",
                                                       "0 < score - Nanopore score < 0.5",
                                                       "score - Nanopore score >= 0.5")))
png(gsub("\\.rds$", "_summary.png", outrds), width = 8, height = 8,
    unit = "in", res = 300)
ggplot(scorecomp, aes(x = method, y = nbr_genes, fill = category)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("blue", "green", "red", "cyan", "pink"), name = "") + 
  theme_bw() + xlab("") + ylab("Number or genes") + 
  guides(fill = guide_legend(nrow = 3)) + theme(legend.position = "bottom")
dev.off()

## Write genes with lower score in nanopore to files
for (cn in setdiff(colnames(scoresub), targetmethod)) {
  tmp <- scoresub[order(scoresub[, targetmethod] - scoresub[, cn]), ]
  genes_to_write <- rownames(tmp)[which(tmp[, targetmethod] < tmp[, cn])]
  write.table(scores %>% dplyr::filter(gene %in% genes_to_write & 
                                         method %in% c(cn, targetmethod)) %>%
                dplyr::select(gene, method, count, length_diff_3putrs_samestart,
                              exoncount, introncount, intron_exon_ratio,
                              uniqjuncreads, mmjuncreads, uniqjuncfraction,
                              nbr_junctions_in_gene, score) %>%
                tidyr::gather(variable, value, count, score) %>%
                tidyr::unite(temp, method, variable) %>%
                tidyr::spread(temp, value) %>%
                dplyr::mutate(gene = factor(gene, levels = genes_to_write)) %>%
                dplyr::arrange(gene),
              file = gsub("\\.rds$", paste0("_", targetmethod, "_lowerthan_", cn, ".txt"), outrds),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

saveRDS(NULL, file = outrds)
date()
sessionInfo()
