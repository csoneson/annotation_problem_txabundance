args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

quantmethods <- strsplit(quantmethods, ",")[[1]]

print(combcovrds)
print(quantmethods)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## Read combined coverage table
x <- readRDS(combcovrds)
x <- x$jcovscaled

## Extract required information
y <- x %>% dplyr::select(gene, junctionid, method, pred.cov, scaled.cov, 
                         uniqreads, mmreads) %>%
  tidyr::gather(covtype, coverage, pred.cov, scaled.cov) %>%
  dplyr::mutate(fracmm = mmreads/(uniqreads+mmreads)) %>%
  dplyr::mutate(fracmm = replace(fracmm, is.na(fracmm), 0)) %>%
  dplyr::mutate(covtype = replace(covtype, covtype == "pred.cov", "Predicted coverage"),
                covtype = replace(covtype, covtype == "scaled.cov", "Scaled predicted coverage")) %>%
  dplyr::filter(method %in% quantmethods)

## Calculate correlations
corrs <- y %>% dplyr::group_by(method, covtype) %>% 
  dplyr::summarize(pearson = signif(cor(uniqreads, coverage, method = "pearson", 
                                        use = "pairwise.complete.obs"), 3),
                   spearman = signif(cor(uniqreads, coverage, method = "spearman", 
                                         use = "pairwise.complete.obs"), 3))

png(gsub("rds$", "png", outrds), width = 7, height = 15, unit = "in", res = 300)
ggplot(y, aes(x = uniqreads, y = coverage)) + 
  geom_abline(intercept = 0, slope = 1, color = "black") + 
  geom_point(alpha = 0.3, size = 0.3, aes(color = (fracmm > 0.5))) + 
  geom_label(data = corrs, x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1, 
             aes(label = paste0("Pearson: ", pearson, "\nSpearman: ", spearman))) + 
  facet_grid(method ~ covtype) + 
  xlab("Number of uniquely mapping reads spanning junction") + 
  ylab("Predicted number of reads spanning junction") + 
  scale_color_manual(name = "Fraction\nmultimapping\nreads > 0.5", 
                     values = c(`TRUE` = "red", `FALSE` = "blue")) + 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) + 
  theme_bw()
dev.off()

saveRDS(NULL, outrds)

date()
sessionInfo()
