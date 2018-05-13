args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

quantmethods <- strsplit(quantmethods, ",")[[1]]

print(scorerds)
print(quantmethods)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(GGally)
})

## Define new panel function for ggpairs, to show both Pearson and Spearman correlation
combinecor <- function(data, mapping, color = I("black"), sizeRange = c(1, 5), ...) {
  
  # get the x and y data to use the other code
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  
  ct1 <- cor(x, y, method = "pearson", use = "pairwise.complete.obs")
  ct2 <- cor(x, y, method = "spearman", use = "pairwise.complete.obs")
  
  # since we can't print it to get the strsize, just use the max size range
  cex <- max(sizeRange)
  
  
  # plot the cor value
  ggally_text(
    label = paste0("Pearson: \n", signif(ct1, 3), "\nSpearman: \n", signif(ct2, 3)), 
    mapping = mapping,
    xP = 0.5, yP = 0.5, 
    color = color,
    ...
  ) 
}

## Read scores
scores <- readRDS(scorerds)$genes

## Filter
## For consistency, expression filtering is done using Salmon counts
highexpression <- scores %>% dplyr::filter(method == "Salmon" & count > 1000) %>% 
  dplyr::pull(gene)
scores_wide <- scores %>% dplyr::filter(method %in% quantmethods & uniqjuncfraction >= 0.75 & 
                                          gene %in% highexpression) %>%
  dplyr::select(gene, method, score) %>% 
  tidyr::spread(method, score) %>% as.data.frame() %>% tibble::column_to_rownames("gene")

## Pairs plot
png(gsub("\\.rds$", "_pairs.png", outrds), width = 8, height = 8, unit = "in", res = 300)
print(ggpairs(scores_wide,
              lower = list(continuous = wrap("points", alpha = 0.3, size = 0.25)),
              upper = list(continuous = combinecor)) + 
        theme_bw())
dev.off()

## Dendrogram
png(gsub("\\.rds$", "_clustering.png", outrds), width = 5, height = 5, unit = "in", res = 300)
cormat <- cor(scores_wide, use = "pairwise.complete.obs")
distmat <- as.dist(sqrt(2*(1 - cormat)))
hcl <- hclust(distmat)
nodePar <- list(pch = c(NA, 19), col = "black")
plot(as.dendrogram(hcl), type = "rectangle", ylab = "Dissimilarity", nodePar = nodePar)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
