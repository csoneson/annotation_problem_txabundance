args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(combcovrds1)
print(combcovrds2)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

combcov1 <- readRDS(combcovrds1)$jcovscaled
combcov2 <- readRDS(combcovrds2)$jcovscaled

combcov1 <- combcov1 %>% dplyr::mutate(diffcov1 = (scaledcoverage - uniqreads)/(uniqreads + 1)) %>%
  dplyr::rename(score1 = score) %>% 
  dplyr::select(seqnames, start, end, strand, gene, method, diffcov1, score1)
combcov2 <- combcov2 %>% dplyr::mutate(diffcov2 = (scaledcoverage - uniqreads)/(uniqreads + 1)) %>%
  dplyr::rename(score2 = score) %>% 
  dplyr::select(seqnames, start, end, strand, gene, method, diffcov2, score2)

combcov <- dplyr::inner_join(combcov1, combcov2, 
                             by = c("seqnames", "start", "end", "strand", "gene", "method"))

png(gsub("\\.rds$", "_score.png", outrds), width = 8, height = 8,
    units = "in", res = 300)

print(ggplot(combcov %>% dplyr::select(gene, method, score1, score2) %>% distinct(), 
             aes(x = score1, y = score2)) + 
        geom_point(size = 0.5) + facet_wrap(~ method) + theme_bw() + 
        geom_abline(slope = 1, intercept = 0) + xlab("Gene score, sample 1") + 
        ylab("Gene score, sample 2"))

dev.off()

png(gsub("\\.rds$", "_deviations.png", outrds), width = 8, height = 8,
    units = "in", res = 300)

print(ggplot(combcov, aes(x = pmax(abs(diffcov1), abs(diffcov2)), 
                          y = diffcov1 - diffcov2)) + 
        geom_point(size = 0.5) + facet_wrap(~ method) + theme_bw() + 
        geom_abline(slope = 0, intercept = 0) + xlab("Maximal deviation") + 
        ylab("Difference between deviations"))

dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
