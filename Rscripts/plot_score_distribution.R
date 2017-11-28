args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(genesummaryrds)  ## gene summary information
print(genemodels)  ## gene models etc (output of alpine_prepare_for_comparison.R)
print(outrds)  ## output file

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

gsr <- readRDS(genesummaryrds)
x <- readRDS(genemodels)
jcovscaled <- x$jcovscaled %>% dplyr::filter(!is.na(gene))

## Create gene information table
gene_groups <- gsr %>% dplyr::mutate(all_genes = 1) %>% 
  dplyr::mutate(high_expression = as.numeric(count > 1000)) %>%
  dplyr::mutate(length_diff_utr = as.numeric(gene %in% genes_with_lengthdiff))

## Get the score for each gene
gene_scores <- jcovscaled %>% dplyr::select(gene, method, score) %>% 
  dplyr::distinct()

## Put information together
gene_plot <- gene_scores %>% dplyr::left_join(
  gene_groups %>% dplyr::select(gene, all_genes, high_expression, length_diff_utr) %>%
    tidyr::gather(group, included, -gene) %>% dplyr::filter(included == 1)
)

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
             aes(x = method, y = score, color = method)) + geom_violin() + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
        ggtitle("Score distribution") + facet_wrap(~ group) +  
        scale_color_manual(values = c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                                      "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                                      "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                                      "#90C987")[seq_len(length(unique(jcovscaled$method)))]))
print(ggplot(gene_plot, 
             aes(x = group, y = score, color = method)) + geom_violin() + theme_bw() + 
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
