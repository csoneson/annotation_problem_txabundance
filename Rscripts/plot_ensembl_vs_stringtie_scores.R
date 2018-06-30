################################################################################
##                                                                            ##
## Compare scores obtained with different annotations as basis                ##
##                                                                            ##
## Inputs:                                                                    ##
## * scorerdsensembl: list containing abundance estimates and characteristics ##
##                    for junctions, transcripts and genes, as well as gene   ##
##                    scores, for the Ensembl annotation.                     ##
## * scorerdsstringtie: list containing abundance estimates and               ##
##                      characteristics for junctions, transcripts and genes, ##
##                      as well as gene scores, for the StringTie annotation. ##
## * convtablestringtie: conversion table for StringTie genes.                ##
## * convtablestringtietx: conversion table for StringTie transcripts.        ##
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

print(scorerdsensembl)
print(scorerdsstringtie)
print(convtablestringtie)
print(convtablestringtietx)
print(uniqjuncreadsthreshold)
print(uniqjuncfracthreshold)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

stringtie <- readRDS(scorerdsstringtie)
ensembl <- readRDS(scorerdsensembl)
stringtieconversion <- readRDS(convtablestringtie)
stringtieconversiontx <- readRDS(convtablestringtietx)

## Keep only StringTie transcripts with an Ensembl correspondence
stringtieconversiontx <- subset(stringtieconversiontx, !is.na(symbol))

## Plot distribution of all scores
combined <- rbind(stringtie$genes %>% dplyr::select(gene, method, score) %>% 
                    dplyr::mutate(annotation = "StringTie"),
                  ensembl$genes %>% dplyr::select(gene, method, score) %>% 
                    dplyr::filter(!(method %in% c("SalmonCDS", "SalmonKeepDup"))) %>%
                    dplyr::mutate(annotation = "Ensembl"))

png(gsub("\\.rds$", "_distribution_allgenes.png", outrds), width = 8, height = 6, 
    unit = "in", res = 400)
ggplot(combined, aes(x = score, color = annotation)) + 
  geom_line(stat = "density", size = 1.5) + facet_wrap(~ method) + theme_bw() + 
  scale_color_manual(values = c(StringTie = "#B17BA6", Ensembl = "#90C987"))
dev.off()

## Same, but only for genes with enough (unique) junction reads
combinedfilt <- rbind(stringtie$genes %>% 
                        dplyr::filter(uniqjuncreads >= uniqjuncreadsthreshold & 
                                        uniqjuncfraction >= uniqjuncfracthreshold) %>% 
                        dplyr::select(gene, method, score) %>% 
                        dplyr::mutate(annotation = "StringTie"),
                      ensembl$genes %>% 
                        dplyr::filter(uniqjuncreads >= uniqjuncreadsthreshold & 
                                        uniqjuncfraction >= uniqjuncfracthreshold) %>% 
                        dplyr::select(gene, method, score) %>% 
                        dplyr::filter(!(method %in% c("SalmonCDS", "SalmonKeepDup"))) %>%
                        dplyr::mutate(annotation = "Ensembl"))

ggde <- ggplot(combinedfilt, aes(x = score, color = annotation)) + 
  geom_line(stat = "density", size = 1.5) + facet_wrap(~ method) + theme_bw() + 
  scale_color_manual(values = c(StringTie = "#B17BA6", Ensembl = "#90C987")) + 
  xlab("Score") + ylab("Density") + theme(legend.position = "bottom")

## Keep only genes with enough (unique) junction reads with both annotations,
## and compare the scores of the genes that can be matched
a <- stringtie$genes %>% 
  dplyr::filter(uniqjuncreads >= uniqjuncreadsthreshold & 
                  uniqjuncfraction >= uniqjuncfracthreshold) %>% 
  dplyr::select(gene, method, score) %>%
  dplyr::left_join(stringtieconversion %>% dplyr::select(gene, symbol) %>%
                     dplyr::distinct(), by = "gene") %>%
  dplyr::select(-gene) %>% dplyr::rename(gene = symbol) %>%
  dplyr::rename(stringtie = score) %>%
  dplyr::inner_join(ensembl$genes %>% 
                      dplyr::filter(uniqjuncreads >= uniqjuncreadsthreshold & 
                                      uniqjuncfraction > uniqjuncfracthreshold) %>%
                      dplyr::select(gene, method, score) %>%
                      dplyr::filter(!(method %in% c("SalmonCDS", "SalmonKeepDup"))) %>%
                      dplyr::rename(ensembl = score), by = c("gene", "method")) %>%
  dplyr::mutate(method = factor(method, levels = unique(method)))

ggcor <- ggplot(a, aes(x = stringtie, y = ensembl)) + geom_point(alpha = 0.3) +
  facet_wrap(~ method) + geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + xlab("StringTie score") + ylab("Ensembl score")

ggbw <- ggplot(a %>% dplyr::group_by(method) %>% 
                 dplyr::summarize(stringtieworse = sum(stringtie > ensembl), 
                                  equal = sum(stringtie == ensembl), 
                                  ensemblworse = sum(ensembl > stringtie)) %>% 
                 tidyr::gather(comparison, value, -method) %>%
                  dplyr::mutate(comparison = replace(comparison, comparison == "stringtieworse", "StringTie score > Ensembl score"),
                                comparison = replace(comparison, comparison == "equal", "StringTie score = Ensembl score"),
                                comparison = replace(comparison, comparison == "ensemblworse", "StringTie score < Ensembl score")), 
               aes(x = method, y = value, fill = comparison)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#B17BA6", "#90C987", "#7BAFDE"), name = "", na.value = "grey50") + 
  xlab("") + ylab("Number of genes") + theme(legend.direction = "horizontal",
                                             legend.justification = "center",
                                             legend.box.just = "bottom",
                                             legend.position = "bottom") + 
  guides(fill = guide_legend(nrow = 2))

png(gsub("\\.rds$", "_combined_expressed.png", outrds), width = 10, height = 10,
    unit = "in", res = 400)
plot_grid(ggcor,
          plot_grid(ggde, ggbw, nrow = 1, rel_widths = c(1, 1), labels = c("B", "C")),
          ncol = 1, rel_heights = c(1, 1), labels = c("A", ""))
dev.off()

## For each (shared) gene and method, what fraction of the expression comes from
## "new" StringTie transcripts? And conversely, what fraction of the expression
## from the Ensembl catalog comes from transcripts that are not in the StringTie
## catalog?
strexpr <- stringtie$transcripts %>% dplyr::filter(symbol %in% ensembl$transcripts$gene) %>% 
  dplyr::group_by(symbol, method) %>%
  dplyr::summarize(fracExprSTR = sum(TPM[!(transcript %in% stringtieconversiontx$tx)])/sum(TPM)) %>%
  dplyr::mutate(fracExprSTR = replace(fracExprSTR, is.na(fracExprSTR), 0))
ensexpr <- ensembl$transcripts %>% dplyr::filter(gene %in% stringtie$transcripts$symbol) %>%
  dplyr::group_by(gene, method) %>%
  dplyr::summarize(fracExprRemoved = sum(TPM[!(transcript %in% stringtieconversiontx$symbol)])/sum(TPM)) %>%
  dplyr::mutate(fracExprRemoved = replace(fracExprRemoved, is.na(fracExprRemoved), 0))

b <- a %>% dplyr::left_join(strexpr %>% dplyr::rename(gene = symbol), 
                            by = c("gene", "method")) %>%
  dplyr::left_join(ensexpr, by = c("gene", "method")) %>%
  dplyr::mutate(str_vs_ens = NA_character_) %>%
  dplyr::mutate(str_vs_ens = replace(str_vs_ens, stringtie > ensembl, 
                                     "StringTie score > Ensembl score"),
                str_vs_ens = replace(str_vs_ens, stringtie < ensembl, 
                                     "StringTie score < Ensembl score"),
                str_vs_ens = replace(str_vs_ens, stringtie == ensembl, 
                                     "StringTie score = Ensembl score")) %>%
  dplyr::left_join(stringtie$transcripts %>% dplyr::select(gene, symbol) %>% 
                     dplyr::distinct() %>% 
                     dplyr::rename(strid = gene), by = c("gene" = "symbol"))
bb <- tidyr::gather(b, fractype, fracexpr, fracExprSTR, fracExprRemoved)

method_colors <- c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                   "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                   "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                   "#90C987")[seq_len(length(unique(bb$method)))]
names(method_colors) <- unique(bb$method)

png(gsub("\\.rds$", "_new_missing_tx.png", outrds), width = 12, height = 7,
    unit = "in", res = 400)
ggplot(bb %>% dplyr::mutate(fractype = replace(fractype, fractype == "fracExprRemoved", 
                                               "Removed Ensembl transcripts"), 
                            fractype = replace(fractype, fractype == "fracExprSTR",
                                               "New StringTie transcripts")), 
       aes(x = str_vs_ens, y = fracexpr, color = method)) + 
  geom_boxplot() + 
  facet_wrap(~ fractype) + theme_bw() + xlab("") + 
  theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1)) + 
  ylab("Fraction of gene TPM") + scale_color_manual(values = method_colors, name = "")
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
