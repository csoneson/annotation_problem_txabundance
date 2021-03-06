################################################################################
##                                                                            ##
## Compare scores obtained with different annotations as basis                ##
##                                                                            ##
## Inputs:                                                                    ##
## * scorerdsensembl: list containing abundance estimates and characteristics ##
##                    for junctions, transcripts and genes, as well as gene   ##
##                    scores, for the Ensembl annotation.                     ##
## * scorerdschess: list containing abundance estimates and characteristics   ##
##                  for junctions, transcripts and genes, as well as gene     ##
##                  scores, for the CHESS annotation.                         ##
## * convtablechess: conversion table for CHESS genes.                        ##
## * quantmethods: string containing the quantification methods to consider,  ##
##                 separated by commas (no spaces)                            ##
## * uniqjuncreadsthreshold: the total number of uniquely mapping junction    ##
##                           reads (in a gene), only genes with more than     ##
##                           this number for both annotations will be used    ##
##                           for the comparison                               ##
## * uniqjuncfracthreshold: the fraction of uniquely mapping junction reads.  ##
##                          Only genes with more than this fraction for both  ##
##                          annotations will be used for the comparison       ##
## * outrds: output rds file. The name will be used to determine the name of  ##
##           the output figures.                                              ##
##                                                                            ##
## Outputs:                                                                   ##
## * png figures comparing the scores obtained with the two annotations       ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

quantmethods <- strsplit(quantmethods, ",")[[1]]

print(scorerdsensembl)
print(scorerdschess)
print(convtablechess)
print(quantmethods)
print(uniqjuncreadsthreshold)
print(uniqjuncfracthreshold)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

chess <- readRDS(scorerdschess)
ensembl <- readRDS(scorerdsensembl)
chessconversion <- readRDS(convtablechess)

## Plot distribution of all scores
combined <- rbind(chess$genes %>% dplyr::select(gene, method, score) %>% 
                    dplyr::mutate(annotation = "CHESS"),
                  ensembl$genes %>% dplyr::select(gene, method, score) %>% 
                    dplyr::filter(method %in% c("Salmon", "kallisto")) %>%
                    dplyr::mutate(annotation = "Ensembl")) %>%
  dplyr::filter(method %in% quantmethods)

png(gsub("\\.rds$", "_distribution_allgenes.png", outrds), width = 8, height = 6, 
    unit = "in", res = 400)
ggplot(combined, aes(x = score, color = annotation)) + 
  geom_line(stat = "density", size = 1.5) + facet_wrap(~ method) + theme_bw() + 
  scale_color_manual(values = c(CHESS = "#B17BA6", Ensembl = "#90C987"))
dev.off()

## Same, but only for genes with enough (unique) junction reads
combinedfilt <- rbind(chess$genes %>% 
                        dplyr::filter(uniqjuncreads >= uniqjuncreadsthreshold & 
                                        uniqjuncfraction >= uniqjuncfracthreshold) %>% 
                        dplyr::select(gene, method, score) %>% 
                        dplyr::mutate(annotation = "CHESS"),
                      ensembl$genes %>% 
                        dplyr::filter(uniqjuncreads >= uniqjuncreadsthreshold & 
                                        uniqjuncfraction >= uniqjuncfracthreshold) %>% 
                        dplyr::select(gene, method, score) %>% 
                        dplyr::filter(method %in% c("Salmon", "kallisto")) %>%
                        dplyr::mutate(annotation = "Ensembl")) %>%
  dplyr::filter(method %in% quantmethods)

ggde <- ggplot(combinedfilt, aes(x = score, color = annotation)) + 
  geom_line(stat = "density", size = 1.5) + facet_wrap(~ method) + theme_bw() + 
  scale_color_manual(values = c(CHESS = "#B17BA6", Ensembl = "#90C987"), name = "") + 
  xlab("JCC score") + ylab("Density") + theme(legend.position = "bottom",
                                              axis.title = element_text(size = 15))

## Keep only genes with enough (unique) junction reads with both annotations,
## and compare the scores of the genes that can be matched
a <- chess$genes %>% 
  dplyr::filter(uniqjuncreads >= uniqjuncreadsthreshold & 
                  uniqjuncfraction >= uniqjuncfracthreshold) %>% 
  dplyr::select(gene, method, score) %>%
  dplyr::left_join(chessconversion %>% dplyr::select(gene, symbol) %>%
                     dplyr::distinct(), by = "gene") %>%
  dplyr::select(-gene) %>% dplyr::rename(gene = symbol) %>%
  dplyr::rename(chess = score) %>%
  dplyr::inner_join(ensembl$genes %>% 
                      dplyr::filter(uniqjuncreads >= uniqjuncreadsthreshold & 
                                      uniqjuncfraction > uniqjuncfracthreshold) %>%
                      dplyr::select(gene, method, score) %>%
                      dplyr::filter(method %in% c("Salmon", "kallisto")) %>%
                      dplyr::rename(ensembl = score), by = c("gene", "method")) %>%
  dplyr::filter(method %in% quantmethods) %>%
  dplyr::mutate(method = factor(method, levels = unique(method)))
  
ggcor <- ggplot(a, aes(x = chess, y = ensembl)) + geom_point(alpha = 0.3) +
  facet_wrap(~ method) + geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + xlab("JCC score (CHESS)") + ylab("JCC score (Ensembl)") + 
  theme(axis.title = element_text(size = 15))

ggbw <- ggplot(a %>% dplyr::group_by(method) %>% 
                 dplyr::summarize(chessmuchworse = sum(chess - ensembl > 0.1),
                                  chesslittleworse = sum(chess - ensembl > 0 & 
                                                           chess - ensembl <= 0.1),
                                  equal = sum(chess == ensembl), 
                                  ensembllittleworse = sum(ensembl - chess > 0 & 
                                                             ensembl - chess <= 0.1),
                                  ensemblmuchworse = sum(ensembl - chess > 0.1)) %>% 
                 tidyr::gather(comparison, value, -method) %>%
                 dplyr::mutate(comparison = replace(comparison, comparison == "chessmuchworse", "CHESS score-Ensembl score > 0.1"),
                               comparison = replace(comparison, comparison == "chesslittleworse", "0 < CHESS score-Ensembl score <= 0.1"),
                               comparison = replace(comparison, comparison == "equal", "CHESS score = Ensembl score"),
                               comparison = replace(comparison, comparison == "ensembllittleworse", "-0.1 <= CHESS score-Ensembl score < 0"),
                               comparison = replace(comparison, comparison == "ensemblmuchworse", "CHESS score-Ensembl score < -0.1")) %>%
                 dplyr::mutate(comparison = factor(comparison,
                                                   levels = c("CHESS score-Ensembl score > 0.1",
                                                              "0 < CHESS score-Ensembl score <= 0.1",
                                                              "CHESS score = Ensembl score",
                                                              "-0.1 <= CHESS score-Ensembl score < 0",
                                                              "CHESS score-Ensembl score < -0.1"))), 
               aes(x = method, y = value, fill = comparison)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#1965B0", "#7BAFDE", "#777777", "#90C987", "#4EB265"), 
                    name = "", na.value = "grey50") + 
  xlab("") + ylab("Number of genes") + theme(legend.direction = "horizontal",
                                             legend.justification = "center",
                                             legend.box.just = "bottom",
                                             legend.position = "bottom",
                                             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  guides(fill = guide_legend(nrow = 5))

png(gsub("\\.rds$", "_combined_expressed.png", outrds), width = 12, height = 10,
    unit = "in", res = 400)
plot_grid(ggcor,
          plot_grid(ggde, ggbw, nrow = 1, rel_widths = c(1, 1), labels = c("B", "C")),
          ncol = 1, rel_heights = c(1, 1), labels = c("A", ""))
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
