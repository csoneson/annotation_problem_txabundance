################################################################################
##                                                                            ##
## Compare scores obtained with different annotations as basis                ##
##                                                                            ##
## Inputs:                                                                    ##
## * genecharensembl: data frame with gene characteristics, for the Ensembl   ##
##                    annotation.                                             ##
## * genecharchess: data frame with gene characteristics, for the CHESS       ##
##                  annotation.                                               ##
## * convtablechess: conversion table for CHESS genes.                        ##
## * outrds: output rds file. The name will be used to determine the name of  ##
##           the output figures.                                              ##
##                                                                            ##
## Outputs:                                                                   ##
## * png figures comparing characteristics of the two annotation catalogs     ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(genecharensembl)
print(genecharchess)
print(convtablechess)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

chess <- readRDS(genecharchess)
ensembl <- readRDS(genecharensembl)
chessconversion <- readRDS(convtablechess)

df <- rbind(
  chess %>% 
    dplyr::select(gene_id, nbr_transcripts, median_transcript_length,
                  max_transcript_length, min_transcript_length,
                  median_nbr_exons_per_tx, max_nbr_exons_per_tx,
                  min_nbr_exons_per_tx) %>% 
    dplyr::mutate(annotation = "CHESS"),
  ensembl %>% 
    dplyr::select(gene_id, nbr_transcripts, median_transcript_length,
                  max_transcript_length, min_transcript_length,
                  median_nbr_exons_per_tx, max_nbr_exons_per_tx,
                  min_nbr_exons_per_tx) %>% 
    dplyr::mutate(annotation = "Ensembl"))

## Number of transcripts per gene
png(gsub("\\.rds$", "_nbrtxpergene.png", outrds), width = 8, height = 6, 
    unit = "in", res = 400)
ggplot(df %>% dplyr::mutate(nbr_transcripts = replace(nbr_transcripts, nbr_transcripts > 30, 30)), 
       aes(x = nbr_transcripts, fill = annotation)) +
  geom_bar(position = "dodge") + scale_y_sqrt() + 
  xlab("Number of transcripts per gene") + ylab("Number of genes") + 
  scale_fill_manual(values = c(CHESS = "#B17BA6", Ensembl = "#90C987"), name = "") + 
  theme_bw() + theme(legend.position = "bottom")
dev.off()

## Median transcript length
png(gsub("\\.rds$", "_mediantxlength.png", outrds), width = 8, height = 6, 
    unit = "in", res = 400)
ggplot(df, aes(x = median_transcript_length, color = annotation)) +
  geom_line(stat = "density", size = 1.5) + scale_x_sqrt() + 
  xlab("Median transcript length per gene") + ylab("Density") + 
  scale_color_manual(values = c(CHESS = "#B17BA6", Ensembl = "#90C987"), name = "") + 
  theme_bw() + theme(legend.position = "bottom")
dev.off()

## Min transcript length
png(gsub("\\.rds$", "_mintxlength.png", outrds), width = 8, height = 6, 
    unit = "in", res = 400)
ggplot(df, aes(x = min_transcript_length, color = annotation)) +
  geom_line(stat = "density", size = 1.5) + scale_x_sqrt() + 
  xlab("Minimum transcript length per gene") + ylab("Density") + 
  scale_color_manual(values = c(CHESS = "#B17BA6", Ensembl = "#90C987"), name = "") + 
  theme_bw() + theme(legend.position = "bottom")
dev.off()

## Max transcript length
png(gsub("\\.rds$", "_maxtxlength.png", outrds), width = 8, height = 6, 
    unit = "in", res = 400)
ggplot(df, aes(x = max_transcript_length, color = annotation)) +
  geom_line(stat = "density", size = 1.5) + scale_x_sqrt() + 
  xlab("Maximum transcript length per gene") + ylab("Density") + 
  scale_color_manual(values = c(CHESS = "#B17BA6", Ensembl = "#90C987"), name = "") + 
  theme_bw() + theme(legend.position = "bottom")
dev.off()

## Number of transcripts for the "same" gene in the two annotations
df2 <- chess %>% 
  dplyr::left_join(chessconversion %>% dplyr::select(gene, symbol) 
                   %>% dplyr::distinct(), by = c("gene_id" = "gene")) %>%
  dplyr::mutate(nbr_transcripts_chess = nbr_transcripts,
                median_transcript_length_chess = median_transcript_length) %>%
  dplyr::select(symbol, nbr_transcripts_chess, median_transcript_length_chess) %>%
  dplyr::inner_join(ensembl %>% 
                      dplyr::mutate(nbr_transcripts_ensembl = nbr_transcripts,
                                    median_transcript_length_ensembl = median_transcript_length) %>%
                      dplyr::select(gene_id, nbr_transcripts_ensembl, median_transcript_length_ensembl),
                    by = c("symbol" = "gene_id"))
png(gsub("\\.rds$", "_nbrtxpergene_shared_genes.png", outrds),
    width = 8, height = 6, unit = "in", res = 400)
ggplot(df2, aes(x = nbr_transcripts_chess, y = nbr_transcripts_ensembl)) + 
  geom_abline(slope = 1, intercept = 0) + geom_point(alpha = 0.3, size = 1) + 
  theme_bw() + xlab("Number of transcripts per gene (CHESS)") + 
  ylab("Number of transcripts per gene (Ensembl)")
dev.off()

png(gsub("\\.rds$", "_mediantxlength_shared_genes.png", outrds),
    width = 8, height = 6, unit = "in", res = 400)
ggplot(df2, aes(x = median_transcript_length_chess, y = median_transcript_length_ensembl)) + 
  geom_abline(slope = 1, intercept = 0) + geom_point(alpha = 0.3, size = 1) + 
  theme_bw() + xlab("Median transcript length per gene (CHESS)") + 
  ylab("Median transcript length per gene (Ensembl)")
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
