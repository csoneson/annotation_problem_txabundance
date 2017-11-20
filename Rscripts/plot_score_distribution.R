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
jcovscaled <- x$jcovscaled

pdf(gsub("rds$", "pdf", outrds))
tmp1 <- gsr %>% dplyr::filter(maxnbrex > 1 & 
                                count > 1000 & 
                                maxlength > 350)
tmp1 <- jcovscaled %>% dplyr::filter(gene %in% tmp1$gene)
print(ggplot(tmp1 %>% dplyr::select(gene, method, score) %>% dplyr::distinct(), 
             aes(x = score, color = method)) + geom_density() + theme_bw() + 
        ggtitle("Score distribution, max nbr exons > 1, count > 1000, max length > 350") + 
        scale_color_manual(values = c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                                      "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                                      "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                                      "#90C987")[seq_len(length(unique(jcovscaled$method)))]))
print(ggplot(tmp1 %>% dplyr::select(gene, method, score) %>% dplyr::distinct(), 
             aes(x = method, y = score, color = method)) + geom_boxplot() + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
        ggtitle("Score distribution, max nbr exons > 1, count > 1000, max length > 350") + 
        scale_color_manual(values = c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                                      "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                                      "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                                      "#90C987")[seq_len(length(unique(jcovscaled$method)))]))
## Rank distribution (1 is best)
print(ggplot(tmp1 %>% dplyr::select(gene, method, score) %>% dplyr::distinct() %>% 
               dplyr::mutate(score = replace(score, is.na(score), 10)) %>% 
               dplyr::group_by(gene) %>% dplyr::mutate(rank = rank(score)) %>%
               dplyr::mutate(keep = !(var(rank) == 0)) %>% dplyr::filter(keep), 
             aes(x = rank)) + geom_bar() + facet_wrap(~method) + theme_bw() + 
        ggtitle("Rank distribution, max nbr exons > 1, count > 1000, max length > 350"))

dev.off()

saveRDS(NULL, file = outrds)
sessionInfo()
date()
