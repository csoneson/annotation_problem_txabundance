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
## * A png figure comparing the scores of all methods to those of the target  ##
##   method                                                                   ##
## * A png figure summarizing the number of genes for which the target method ##
##   gives higher/lower/equal score compared to each other method             ##
## * For each method, a text file with the genes for which the target method  ##
##   gave a lower score                                                       ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(scorerdsensembl)
print(scorerdschess)
print(convtablechess)
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
                    dplyr::mutate(annotation = "Ensembl"))

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
                        dplyr::mutate(annotation = "Ensembl"))

png(gsub("\\.rds$", "_distribution_expressed.png", outrds), width = 8, height = 6, 
    unit = "in", res = 400)
ggplot(combinedfilt, aes(x = score, color = annotation)) + 
  geom_line(stat = "density", size = 1.5) + facet_wrap(~ method) + theme_bw() + 
  scale_color_manual(values = c(CHESS = "#B17BA6", Ensembl = "#90C987"))
dev.off()

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
  dplyr::mutate(method = factor(method, levels = unique(method)))
  
png(gsub("\\.rds$", "_correlation_expressed.png", outrds), width = 8, height = 6,
    unit = "in", res = 400)
ggplot(a, aes(x = chess, y = ensembl)) + geom_point(alpha = 0.3) +
  facet_wrap(~ method) + geom_abline(slope = 1, intercept = 0) + 
  theme_bw()
dev.off()

png(gsub("\\.rds$", "_betterworse_expressed.png", outrds), width = 8, height = 6,
    unit = "in", res = 400)
ggplot(as.data.frame(table(a$method, chessworse = a$chess > a$ensembl, useNA = "ifany")) %>%
         dplyr::mutate(chessworse = as.character(chessworse)) %>% 
         dplyr::mutate(chessworse = replace(chessworse, chessworse == "TRUE", "CHESS score > Ensembl score"),
                       chessworse = replace(chessworse, chessworse == "FALSE", "CHESS score < Ensembl score")), 
       aes(x = Var1, y = Freq, fill = chessworse)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#B17BA6", "#90C987"), name = "", na.value = "grey50") + 
  xlab("") + ylab("Number of genes")
dev.off()

date()
sessionInfo()
