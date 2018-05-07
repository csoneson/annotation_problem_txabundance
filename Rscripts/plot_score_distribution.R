args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

quantmethods <- strsplit(quantmethods, ",")[[1]]

print(scorerds)
print(gexrds)  ## gene expression
print(geneinfords)  ## gene information
print(exoncountstxt)
print(introncountstxt)
print(quantmethods)
print(outrds)  ## output file

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

scores <- readRDS(scorerds)
gex <- readRDS(gexrds)
geneinfo <- readRDS(geneinfords)
exoncounts <- read.delim(exoncountstxt, skip = 1, header = TRUE, as.is = TRUE) %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
  setNames(c("gene", "exoncount"))
introncounts <- read.delim(introncountstxt, skip = 1, header = TRUE, as.is = TRUE) %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
  setNames(c("gene", "introncount"))

## Define colors
method_colors <- c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                   "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                   "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                   "#90C987")[seq_len(length(unique(scores$method)))]
names(method_colors) <- unique(scores$method)

## Calculate the fraction of multimapping junction reads for each gene
mmfrac <- scores %>% dplyr::select(junctionid, gene, uniqreads, mmreads) %>%
  dplyr::distinct() %>% dplyr::group_by(gene) %>% 
  dplyr::summarize(mmreads = sum(mmreads), uniqreads = sum(uniqreads)) %>% 
  dplyr::mutate(mmfraction = mmreads/(mmreads + uniqreads)) %>%
  dplyr::select(gene, mmfraction, uniqreads, mmreads) %>%
  dplyr::rename(uniq_junc_reads = uniqreads, mm_junc_reads = mmreads)

## Get the score for each gene
gene_scores <- list()
gene_scores[["score"]] <- scores %>% dplyr::select(gene, method, score) %>% 
  dplyr::distinct() %>%
  dplyr::left_join(mmfrac, by = "gene") %>%
  dplyr::left_join(geneinfo, by = c("gene" = "gene_id")) %>%
  dplyr::left_join(gex) %>%
  dplyr::group_by(gene) %>%
  dplyr::left_join(gex %>% dplyr::filter(method == "Salmon") %>% 
                     dplyr::select(gene, count) %>% dplyr::rename(salmon_count = count)) %>% 
  dplyr::filter(method %in% quantmethods)

gene_scores[["scoreMM"]] <- scores %>% dplyr::select(gene, method, scoreMM) %>% 
  dplyr::distinct() %>%
  dplyr::left_join(mmfrac, by = "gene") %>%
  dplyr::left_join(geneinfo, by = c("gene" = "gene_id")) %>%
  dplyr::left_join(gex) %>%
  dplyr::group_by(gene) %>%
  dplyr::left_join(gex %>% dplyr::filter(method == "Salmon") %>% 
                     dplyr::select(gene, count) %>% dplyr::rename(salmon_count = count)) %>% 
  dplyr::filter(method %in% quantmethods)

## Add ratio between intron and exon reads
intron_exon_ratio <- dplyr::full_join(exoncounts, introncounts) %>% 
  dplyr::mutate(introncount = replace(introncount, is.na(introncount), 0)) %>%
  dplyr::mutate(intron_exon_ratio = introncount/exoncount) %>%
  dplyr::mutate(intron_exon_ratio = replace(intron_exon_ratio, exoncount==0 & introncount==0, 0))
gene_scores <- lapply(gene_scores, function(w) {
  dplyr::left_join(w, intron_exon_ratio)
})

for (sc in c("score", "scoreMM")) {
  png(gsub("\\.rds$", paste0("_", sc, ".png"), outrds), width = 14, 
      height = 14, unit = "in", res = 300)
  
  ## Plot score distribution
  print(
    cowplot::plot_grid(
      ggplot(gene_scores[[sc]], aes_string(x = 1, y = sc, fill = "method")) + 
        geom_violin() + theme_bw() + facet_wrap(~ method, nrow = 1) + 
        scale_color_manual(values = method_colors) + 
        theme(legend.position = "none") + xlab(""), 
      
      ggplot(gene_scores[[sc]], aes_string(x = "mmfraction", y = sc, color = "method")) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) + 
        geom_smooth(color = "black", method = "loess") + xlab("Fraction multimapping reads") + 
        theme(legend.position = "none"),
      
      ggplot(gene_scores[[sc]] %>% dplyr::mutate(intron_exon_ratio = replace(intron_exon_ratio, 
                                                                             intron_exon_ratio > 10, 10)), 
             aes_string(x = "intron_exon_ratio", y = sc, color = "method")) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) + 
        geom_smooth(color = "black", method = "loess") + xlab("Intron/exon count ratio") + 
        theme(legend.position = "none"),
      
      ggplot(gene_scores[[sc]], aes_string(x = "salmon_count", y = sc, color = "method")) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) + 
        geom_smooth(color = "black", method = "loess") + xlab("Salmon gene count") + scale_x_sqrt() + 
        theme(legend.position = "none"),
      
      ggplot(gene_scores[[sc]], aes_string(x = "length_diff_3putrs_samestart", y = sc, color = "method")) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) + 
        geom_smooth(color = "black", method = "loess") + xlab("Length difference of 3'UTRs with same start") + 
        theme(legend.position = "none"),
      
      ncol = 1      
    )
  )
  
  dev.off()
}


saveRDS(NULL, file = outrds)
sessionInfo()
date()
