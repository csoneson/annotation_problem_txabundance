args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(combcovrds)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(GGally)
})

jcov <- readRDS(combcovrds)$jcovscaled
combcovall <- jcov %>%
  dplyr::select(gene, method, score) %>% 
  dplyr::distinct() %>%
  tidyr::spread(method, score) %>% 
  as.data.frame()
rownames(combcovall) <- combcovall$gene
combcovall$gene <- NULL

png(gsub("\\.rds$", "_score_correlation_all.png", outrds), width = 8, height = 8,
    units = "in", res = 300)

print(ggpairs(combcovall, 
              lower = list(continuous = wrap("points", alpha = 0.3, size = 0.25))) + 
        theme_bw())

dev.off()

combcovlowmm <- jcov %>%
  dplyr::select(gene, method, uniqreads, mmreads, score) %>%
  dplyr::group_by(gene, method, score) %>% 
  dplyr::summarise(uniqreads = sum(uniqreads), mmreads = sum(mmreads)) %>%
  dplyr::filter(uniqreads/(uniqreads + mmreads) > 0.8) %>%
  dplyr::select(-uniqreads, -mmreads) %>%
  dplyr::distinct() %>%
  tidyr::spread(method, score) %>% 
  as.data.frame()
rownames(combcovlowmm) <- combcovlowmm$gene
combcovlowmm$gene <- NULL

png(gsub("\\.rds$", "_score_correlation_lowmm.png", outrds), width = 8, height = 8,
    units = "in", res = 300)

print(ggpairs(combcovlowmm, 
              lower = list(continuous = wrap("points", alpha = 0.3, size = 0.25))) + 
        theme_bw())

dev.off()

dm <- as.dist(1 - cor(combcovall, use = "pairwise.complete.obs"))
hc <- hclust(dm)

png(gsub("\\.rds$", "_score_correlation_dendrogram.png", outrds), width = 8, height = 8,
    units = "in", res = 300)

plot(hc)

dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
