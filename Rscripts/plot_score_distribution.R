args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(covrds)
print(gexrds)  ## gene expression
print(geneinfords)  ## gene information
print(exoncountstxt)
print(introncountstxt)
print(outrds)  ## output file

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

cov <- readRDS(covrds)
gex <- readRDS(gexrds)
geneinfo <- readRDS(geneinfords)
exoncounts <- read.delim(exoncountstxt, skip = 1, header = TRUE, as.is = TRUE) %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
  setNames(c("gene", "exoncount"))
introncounts <- read.delim(introncountstxt, skip = 1, header = TRUE, as.is = TRUE) %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
  setNames(c("gene", "introncount"))

jcovscaled <- cov$jcovscaled

## Define colors
method_colors <- c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                   "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                   "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                   "#90C987")[seq_len(length(unique(jcovscaled$method)))]
names(method_colors) <- unique(jcovscaled$method)

## Calculate the fraction of multimapping junction reads for each gene
mmfrac <- jcovscaled %>% dplyr::select(junctionid, gene, uniqreads, mmreads) %>%
  dplyr::distinct() %>% dplyr::group_by(gene) %>% 
  dplyr::summarize(mmreads = sum(mmreads), uniqreads = sum(uniqreads)) %>% 
  dplyr::mutate(mmfraction = mmreads/(mmreads + uniqreads)) %>%
  dplyr::select(gene, mmfraction, uniqreads, mmreads) %>%
  dplyr::rename(uniq_junc_reads = uniqreads, mm_junc_reads = mmreads)

## Get the score for each gene
gene_scores <- jcovscaled %>% dplyr::select(gene, method, score) %>% 
  dplyr::distinct() %>%
  dplyr::left_join(mmfrac, by = "gene") %>%
  dplyr::left_join(geneinfo, by = c("gene" = "gene_id")) %>%
  dplyr::left_join(gex) %>%
  dplyr::group_by(gene) %>%
  dplyr::left_join(gex %>% dplyr::filter(method == "Salmon") %>% 
                     dplyr::select(gene, count) %>% dplyr::rename(salmon_count = count))

## Add ratio between intron and exon reads
intron_exon_ratio <- dplyr::full_join(exoncounts, introncounts) %>% 
  dplyr::mutate(introncount = replace(introncount, is.na(introncount), 0)) %>%
  dplyr::mutate(intron_exon_ratio = introncount/exoncount) %>%
  dplyr::mutate(intron_exon_ratio = replace(intron_exon_ratio, exoncount==0 & introncount==0, 0))
gene_scores <- dplyr::left_join(gene_scores, intron_exon_ratio)

pdf(gsub("rds$", "pdf", outrds), width = 10)

## Plot score distribution
print(ggplot(gene_scores, aes(x = 1, y = score, color = method)) + 
        geom_violin() + theme_bw() + facet_wrap(~ method) + 
        scale_color_manual(values = method_colors))

## Plot score distribution by fraction multimapping reads
print(ggplot(gene_scores, aes(x = mmfraction, y = score, color = method)) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method) + 
        geom_smooth(color = "black") + xlab("Fraction multimapping reads") + 
        scale_color_manual(values = method_colors))

print(ggplot(gene_scores %>% dplyr::filter(!is.na(mmfraction)), 
             aes(x = mmfraction > 0.5, y = score, color = method)) + 
        geom_violin() + theme_bw() + facet_wrap(~ method) + 
        xlab("Fraction multimapping reads > 50%") + 
        scale_color_manual(values = method_colors))

## For reads with few multimapping reads, plot score distribution by intron/exon ratio
print(ggplot(gene_scores %>% dplyr::filter(!is.na(mmfraction) & mmfraction <= 0.5 & 
                                             !is.na(intron_exon_ratio)),
             aes(x = intron_exon_ratio, y = score, color = method)) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method) + 
        geom_smooth(color = "black") + xlab("Intron/exon count ratio") + 
        scale_color_manual(values = method_colors))

print(ggplot(gene_scores %>% dplyr::filter(!is.na(mmfraction) & mmfraction <= 0.5 & 
                                             !is.na(intron_exon_ratio)),
             aes(x = intron_exon_ratio > 1, y = score, color = method)) + 
        geom_violin() + theme_bw() + facet_wrap(~ method) + 
        xlab("Intron/exon count ratio > 1") + 
        scale_color_manual(values = method_colors))

## For reads with few multimapping reads and low intron/exon ratio, plot score
## distribution by expression
print(ggplot(gene_scores %>% dplyr::filter(!is.na(mmfraction) & mmfraction <= 0.5 & 
                                             !is.na(intron_exon_ratio) & intron_exon_ratio < 1 & 
                                             !is.na(salmon_count)),
             aes(x = salmon_count, y = score, color = method)) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method) + 
        geom_smooth(color = "black") + xlab("Salmon gene count") + scale_x_sqrt() + 
        scale_color_manual(values = method_colors))

print(ggplot(gene_scores %>% dplyr::filter(!is.na(mmfraction) & mmfraction <= 0.5 & 
                                             !is.na(intron_exon_ratio) & intron_exon_ratio < 1 & 
                                             !is.na(salmon_count)),
             aes(x = salmon_count > 10, y = score, color = method)) + 
        geom_violin() + theme_bw() + facet_wrap(~ method) + 
        xlab("Salmon count > 10") + 
        scale_color_manual(values = method_colors))

## For reads with few multimapping reads, low intron/exon ratio and high
## expression, plot score distribution by length difference of 3' UTRs with same
## start
print(ggplot(gene_scores %>% dplyr::filter(!is.na(mmfraction) & mmfraction <= 0.5 & 
                                             !is.na(intron_exon_ratio) & intron_exon_ratio < 1 & 
                                             !is.na(salmon_count) & salmon_count > 10 & 
                                             !is.na(length_diff_3putrs_samestart)),
             aes(x = length_diff_3putrs_samestart, y = score, color = method)) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method) + 
        geom_smooth(color = "black") + xlab("Length difference of 3'UTRs with same start") + 
        scale_color_manual(values = method_colors))

print(ggplot(gene_scores %>% dplyr::filter(!is.na(mmfraction) & mmfraction <= 0.5 & 
                                             !is.na(intron_exon_ratio) & intron_exon_ratio < 1 & 
                                             !is.na(salmon_count) & salmon_count > 10 & 
                                             !is.na(length_diff_3putrs_samestart)),
             aes(x = length_diff_3putrs_samestart > 1000, y = score, color = method)) + 
        geom_violin() + theme_bw() + facet_wrap(~ method) + 
        xlab("Length difference of 3'UTRs with same start > 1000") + 
        scale_color_manual(values = method_colors))

## For reads with few multimapping reads, low intron/exon ratio and high
## expression, plot score distribution by average number of exons/transcript
print(ggplot(gene_scores %>% dplyr::filter(!is.na(mmfraction) & mmfraction <= 0.5 & 
                                             !is.na(intron_exon_ratio) & intron_exon_ratio < 1 & 
                                             !is.na(salmon_count) & salmon_count > 10 & 
                                             !is.na(ave_nbr_exons)),
             aes(x = ave_nbr_exons, y = score, color = method)) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method) + 
        geom_smooth(color = "black") + xlab("Average number of exons per transcript") + 
        scale_color_manual(values = method_colors))

print(ggplot(gene_scores %>% dplyr::filter(!is.na(mmfraction) & mmfraction <= 0.5 & 
                                             !is.na(intron_exon_ratio) & intron_exon_ratio < 1 & 
                                             !is.na(salmon_count) & salmon_count > 10 & 
                                             !is.na(ave_nbr_exons)),
             aes(x = ave_nbr_exons > 10, y = score, color = method)) + 
        geom_violin() + theme_bw() + facet_wrap(~ method) + 
        xlab("Average number of exons per transcript > 10") + 
        scale_color_manual(values = method_colors))

dev.off()

## Create gene categories
# gene_plot <- gene_scores %>% dplyr::select(gene, method, score, all_genes, high_expression,
#                                            high_expr_length_diff_utr, high_expr_few_multimap, 
#                                            high_expr_many_junctions, many_multimap) %>%
#   tidyr::gather(group, included, -gene, -score, -method) %>% dplyr::filter(included == 1)
# 
# 
# pdf(gsub("rds$", "pdf", outrds), width = 14)
# ## Rank distribution (1 is best)
# print(ggplot(gene_plot %>% 
#                dplyr::mutate(score = replace(score, is.na(score), 10)) %>% 
#                dplyr::group_by(gene, group) %>% dplyr::mutate(rank = rank(score)) %>%
#                dplyr::mutate(keep = !(var(rank) == 0)) %>% dplyr::filter(keep), 
#              aes(x = rank)) + geom_bar() + 
#         facet_grid(group ~ method, scales = "free_y") + theme_bw() + 
#         ggtitle("Rank distribution"))
# 
# dev.off()

saveRDS(NULL, file = outrds)
sessionInfo()
date()
