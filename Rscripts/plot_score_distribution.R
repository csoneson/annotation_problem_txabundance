args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(covrds)
print(gexrds)  ## gene expression
print(geneinfords)  ## gene information
print(outrds)  ## output file

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

cov <- readRDS(covrds)
gex <- readRDS(gexrds)
geneinfo <- readRDS(geneinfords)

jcovscaled <- cov$jcovscaled

## Calculate the fraction of multimapping junction reads for each gene
mmfrac <- jcovscaled %>% dplyr::select(junctionid, gene, uniqreads, mmreads) %>%
  dplyr::distinct() %>% dplyr::group_by(gene) %>% 
  dplyr::summarize(mmreads = sum(mmreads), uniqreads = sum(uniqreads)) %>% 
  dplyr::mutate(mmfraction = mmreads/(mmreads + uniqreads)) %>%
  dplyr::select(gene, mmfraction)

## Get the score for each gene
gene_scores <- jcovscaled %>% dplyr::select(gene, method, score) %>% 
  dplyr::distinct() %>%
  dplyr::left_join(mmfrac, by = "gene") %>%
  dplyr::left_join(geneinfo, by = c("gene" = "gene_id")) %>%
  dplyr::left_join(gex) %>%
  dplyr::group_by(gene) %>%
  dplyr::left_join(gex %>% dplyr::filter(method == "Salmon") %>% 
                     dplyr::select(gene, count) %>% dplyr::rename(salmon_count = count))

## Create gene categories
gene_scores <- gene_scores %>% dplyr::mutate(all_genes = 1) %>%
  dplyr::mutate(many_multimap = as.numeric(mmfraction > 0.5)) %>%
  dplyr::mutate(high_expr_few_multimap = as.numeric(salmon_count > 1000 & 
                                                      mmfraction < 0.5)) %>%
  dplyr::mutate(high_expression = as.numeric(salmon_count > 1000)) %>%
  dplyr::mutate(high_expr_length_diff_utr = as.numeric(salmon_count > 1000 & 
                                                         length_diff_3putrs_samestart > 1000)) %>%
  dplyr::mutate(high_expr_many_junctions = as.numeric(salmon_count > 1000 & ave_nbr_exons > 8))

gene_plot <- gene_scores %>% dplyr::select(gene, method, score, all_genes, high_expression,
                                           high_expr_length_diff_utr, high_expr_few_multimap, 
                                           high_expr_many_junctions, many_multimap) %>%
  tidyr::gather(group, included, -gene, -score, -method) %>% dplyr::filter(included == 1)


pdf(gsub("rds$", "pdf", outrds), width = 14)
print(ggplot(gene_plot, 
             aes(x = score, color = method)) + geom_density() + theme_bw() + 
        ggtitle("Score distribution") + facet_wrap(~ group) + 
        scale_color_manual(values = c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                                      "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                                      "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                                      "#90C987")[seq_len(length(unique(jcovscaled$method)))]))
print(ggplot(gene_plot, 
             aes(x = score, color = group)) + geom_density() + theme_bw() + 
        ggtitle("Score distribution") + facet_wrap(~ method))
print(ggplot(gene_plot, 
             aes(x = method, y = score, color = method)) + geom_boxplot() + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
        ggtitle("Score distribution") + facet_wrap(~ group) +  
        scale_color_manual(values = c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                                      "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                                      "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                                      "#90C987")[seq_len(length(unique(jcovscaled$method)))]))
print(ggplot(gene_plot, 
             aes(x = group, y = score, color = method)) + geom_boxplot() + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
        ggtitle("Score distribution") + facet_wrap(~ method) +  
        scale_color_manual(values = c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                                      "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                                      "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                                      "#90C987")[seq_len(length(unique(jcovscaled$method)))]))
## Rank distribution (1 is best)
print(ggplot(gene_plot %>% 
               dplyr::mutate(score = replace(score, is.na(score), 10)) %>% 
               dplyr::group_by(gene, group) %>% dplyr::mutate(rank = rank(score)) %>%
               dplyr::mutate(keep = !(var(rank) == 0)) %>% dplyr::filter(keep), 
             aes(x = rank)) + geom_bar() + 
        facet_grid(group ~ method, scales = "free_y") + theme_bw() + 
        ggtitle("Rank distribution"))

dev.off()

saveRDS(NULL, file = outrds)
sessionInfo()
date()
