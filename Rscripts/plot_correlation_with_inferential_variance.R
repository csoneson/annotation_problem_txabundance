args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(scorerds)
print(salmondir)
print(tx2gene)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(tximport)
})

## Read scores
scores <- readRDS(scorerds)$genes

## Read bootstrap counts on transcript level
tx2gene <- readRDS(tx2gene)
salmonres <- tximport(paste0(salmondir, "/quant.sf"), type = "salmon", txOut = TRUE)
bootcounts <- salmonres$infReps[[1]]
rownames(bootcounts) <- rownames(salmonres$counts)

## Summarize counts on gene level
bootcounts_gene <- as.data.frame(bootcounts) %>% tibble::rownames_to_column("tx") %>%
  dplyr::left_join(tx2gene %>% dplyr::select(tx, gene), by = "tx") %>%
  dplyr::select(-tx) %>% dplyr::group_by(gene) %>% dplyr::summarize_all(funs(sum)) %>%
  as.data.frame() %>% tibble::column_to_rownames(var = "gene")
rownames(bootcounts_gene) <- gsub("\\.[0-9]+$", "", rownames(bootcounts_gene))

## Calculate coefficient of variation for each transcript and gene
cv_tx <- apply(bootcounts, 1, sd)/apply(bootcounts, 1, mean)
cv_tx[is.na(cv_tx)] <- 0
cv_gene <- apply(bootcounts_gene, 1, sd)/apply(bootcounts_gene, 1, mean)
cv_gene[is.na(cv_gene)] <- 0

## Merge with score information. For transcripts, we use the score for the
## corresponding gene
df_tx <- data.frame(tx = names(cv_tx), CV = cv_tx, stringsAsFactors = FALSE) %>%
  dplyr::left_join(tx2gene %>% dplyr::select(tx, gene)) %>%
  dplyr::mutate(tx = gsub("\\.[0-9]+$", "", tx),
                gene = gsub("\\.[0-9]+$", "", gene)) %>% 
  dplyr::inner_join(scores %>% dplyr::filter(method == "Salmon" & count > 1000 & uniqjuncfraction >= 0.75) %>%
                      dplyr::select(gene, score))

df_gene <- scores %>% dplyr::filter(method == "Salmon" & count > 1000 & uniqjuncfraction >= 0.75) %>% 
  dplyr::select(gene, score) %>%
  dplyr::inner_join(data.frame(gene = names(cv_gene), CV = cv_gene, stringsAsFactors = FALSE))

## Plot
png(gsub("rds$", "png", outrds), width = 6, height = 12, unit = "in", res = 300)
cowplot::plot_grid(
  ggplot(df_tx, aes(x = CV, y = score)) + 
    geom_point(alpha = 0.3, size = 1) + geom_smooth() + 
    theme_bw() + ggtitle("Transcript"),
  
  ggplot(df_gene, aes(x = CV, y = score)) + 
    geom_point(alpha = 0.3, size = 1) + geom_smooth(data = df_gene %>% dplyr::filter(CV < 0.04)) + 
    theme_bw() + ggtitle("Gene"),
  
  ncol = 1
)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
