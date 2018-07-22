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
print(scorerdsstringtie)
print(convtablestringtie)
print(convtablestringtietx)
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

## Load method colors
source("Rscripts/define_plot_colors.R")

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
                    dplyr::mutate(annotation = "Ensembl")) %>%
  dplyr::filter(method %in% quantmethods)

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
                        dplyr::mutate(annotation = "Ensembl")) %>%
  dplyr::filter(method %in% quantmethods)

ggde <- ggplot(combinedfilt, aes(x = score, color = annotation)) + 
  geom_line(stat = "density", size = 1.5) + facet_wrap(~ method) + theme_bw() + 
  scale_color_manual(values = c(StringTie = "#B17BA6", Ensembl = "#90C987"),
                     name = "") + 
  xlab("JCC score") + ylab("Density") + 
  theme(legend.position = "bottom",
        axis.title = element_text(size = 15))

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
  dplyr::filter(method %in% quantmethods) %>%
  dplyr::mutate(method = factor(method, levels = unique(method)))

ggcor <- ggplot(a, aes(x = stringtie, y = ensembl)) + geom_point(alpha = 0.3) +
  facet_wrap(~ method) + geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + xlab("JCC score (StringTie)") + ylab("JCC score (Ensembl)") + 
  theme(axis.title = element_text(size = 15))

ggbw <- ggplot(a %>% dplyr::group_by(method) %>% 
                 dplyr::summarize(stringtiemuchworse = sum(stringtie - ensembl > 0.1),
                                  stringtielittleworse = sum(stringtie - ensembl > 0 &
                                                               stringtie - ensembl <= 0.1),
                                  equal = sum(stringtie == ensembl), 
                                  ensembllittleworse = sum(ensembl - stringtie > 0 & 
                                                             ensembl - stringtie <= 0.1),
                                  ensemblmuchworse = sum(ensembl - stringtie > 0.1)) %>% 
                 tidyr::gather(comparison, value, -method) %>%
                  dplyr::mutate(comparison = replace(comparison, comparison == "stringtiemuchworse", "StringTie score-Ensembl score > 0.1"),
                                comparison = replace(comparison, comparison == "stringtielittleworse", "0 < StringTie score-Ensembl score <= 0.1"),
                                comparison = replace(comparison, comparison == "equal", "StringTie score = Ensembl score"),
                                comparison = replace(comparison, comparison == "ensembllittleworse", "-0.1 <= StringTie score-Ensembl score < 0"),
                                comparison = replace(comparison, comparison == "ensemblmuchworse", "StringTie score-Ensembl score < -0.1")) %>%
                 dplyr::mutate(comparison = factor(comparison,
                                                   levels = c("StringTie score-Ensembl score > 0.1",
                                                              "0 < StringTie score-Ensembl score <= 0.1",
                                                              "StringTie score = Ensembl score",
                                                              "-0.1 <= StringTie score-Ensembl score < 0",
                                                              "StringTie score-Ensembl score < -0.1"))), 
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
          plot_grid(ggde, ggbw, nrow = 1, rel_widths = c(1, 1.2), labels = c("B", "C")),
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
  dplyr::mutate(str_vs_ens = replace(str_vs_ens, stringtie - ensembl > 0.1, 
                                     "StringTie score-Ensembl score > 0.1"),
                str_vs_ens = replace(str_vs_ens, stringtie - ensembl > 0 & 
                                       stringtie - ensembl <= 0.1, 
                                     "0 < StringTie score-Ensembl score <= 0.1"),
                str_vs_ens = replace(str_vs_ens, stringtie - ensembl < -0.1, 
                                     "StringTie score-Ensembl score < -0.1"),
                str_vs_ens = replace(str_vs_ens, stringtie - ensembl >= -0.1 & 
                                       stringtie - ensembl < 0,
                                     "-0.1 <= StringTie score-Ensembl score < 0"),
                str_vs_ens = replace(str_vs_ens, stringtie == ensembl, 
                                     "StringTie score = Ensembl score")) %>%
  dplyr::mutate(str_vs_ens = factor(str_vs_ens, 
                                    levels = c("StringTie score-Ensembl score > 0.1",
                                               "0 < StringTie score-Ensembl score <= 0.1",
                                               "StringTie score = Ensembl score",
                                               "-0.1 <= StringTie score-Ensembl score < 0",
                                               "StringTie score-Ensembl score < -0.1"))) %>%
  dplyr::filter(method %in% quantmethods) %>%
  dplyr::left_join(stringtie$transcripts %>% dplyr::select(gene, symbol) %>% 
                     dplyr::distinct() %>% 
                     dplyr::rename(strid = gene), by = c("gene" = "symbol"))
bb <- tidyr::gather(b, fractype, fracexpr, fracExprSTR, fracExprRemoved)

png(gsub("\\.rds$", "_new_missing_tx.png", outrds), width = 15, height = 7,
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
