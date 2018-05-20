args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

quantmethods <- strsplit(quantmethods, ",")[[1]]

print(scorerds)
print(quantmethods)
print(truthrda)
print(truthmodgenesrds)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(GGally)
})

source("Rscripts/helper_plot_functions.R")

## Read transcript abundance estimates
estimates <- readRDS(scorerds)$transcripts %>% 
  dplyr::select(transcript, gene, count, method) %>%
  dplyr::filter(method %in% quantmethods)

## Read true (simulated) abundance estimates
load(truthrda)  ## counts_matrix
truth <- as.data.frame(counts_matrix) %>% dplyr::select(sample_01) %>%
  tibble::rownames_to_column("transcript") %>% dplyr::rename(truth = sample_01)

## Merge
merged <- estimates %>% tidyr::spread(method, count) %>% 
  dplyr::full_join(truth, by = c("transcript")) %>%
  replace(., is.na(.), 0) %>%
  tidyr::gather(method, count, -transcript, -gene, -truth)

## Add gene annotation for modified transcripts
idx <- grep("utrfrom", merged$transcript)
for (i in idx) {
  merged$gene[i] <- merged$gene[match(strsplit(merged$transcript[i], "_")[[1]][1], 
                                      merged$transcript)]
}

## Add label for genes whose transcripts are modified
truthmodgenes <- readRDS(truthmodgenesrds)
merged <- merged %>% dplyr::mutate(modified_gene = gene %in% truthmodgenes)

## Plot estimated vs true counts
## Transcripts
corrs <- merged %>% dplyr::group_by(method) %>% 
  dplyr::summarize(pearson = signif(cor(truth, count, method = "pearson", 
                                        use = "pairwise.complete.obs"), 3),
                   spearman = signif(cor(truth, count, method = "spearman", 
                                         use = "pairwise.complete.obs"), 3))
png(gsub("\\.rds$", "_transcripts.png", outrds), width = 12, height = 12, 
    unit = "in", res = 300)
print(ggplot(merged, aes(x = truth, y = count)) + 
        geom_abline(slope = 1, intercept = 0) + ggtitle("Transcript") + 
        geom_point(alpha = 0.3, size = 0.5, aes(color = modified_gene)) + theme_bw() + 
        geom_label(data = corrs, x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1, 
                   aes(label = paste0("Pearson: ", pearson, "\nSpearman: ", spearman))) + 
        facet_wrap(~ method) + xlab("True count") + ylab("Estimated count"))
dev.off()

## Genes
mergedgene <- merged %>% dplyr::group_by(method, gene) %>% 
  dplyr::summarize(truth = sum(truth), count = sum(count),
                   modified_gene = unique(modified_gene))
corrsgene <- mergedgene %>% dplyr::group_by(method) %>% 
  dplyr::summarize(pearson = signif(cor(truth, count, method = "pearson", 
                                        use = "pairwise.complete.obs"), 3),
                   spearman = signif(cor(truth, count, method = "spearman", 
                                         use = "pairwise.complete.obs"), 3))
png(gsub("\\.rds$", "_genes.png", outrds), width = 12, height = 12, 
    unit = "in", res = 300)
print(ggplot(mergedgene, 
             aes(x = truth, y = count)) + 
        geom_abline(slope = 1, intercept = 0) + ggtitle("Gene") + 
        geom_point(alpha = 0.3, size = 0.5, aes(color = modified_gene)) + theme_bw() + 
        geom_label(data = corrsgene, x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1, 
                   aes(label = paste0("Pearson: ", pearson, "\nSpearman: ", spearman))) + 
        facet_wrap(~ method) + xlab("True count") + ylab("Estimated count"))
dev.off()

## Pairs plot
png(gsub("\\.rds$", "_transcripts_pairs.png", outrds), width = 12, height = 12,
    unit = "in", res = 300)
toplot <- merged %>% dplyr::select(-gene, -truth) %>% 
  tidyr::spread(method, count) %>% dplyr::select(-transcript) %>%
  dplyr::mutate(modified_gene = as.factor(modified_gene))
print(ggpairs(toplot, 
              columns = which(colnames(toplot) != "modified_gene"), 
              mapping = ggplot2::aes(colour = modified_gene),
              lower = list(continuous = wrap("points", alpha = 0.3, size = 0.25))) + 
        theme_bw())
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()

