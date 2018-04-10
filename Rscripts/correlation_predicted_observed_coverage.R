args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(combcovrds)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

combcov <- readRDS(combcovrds)

png(gsub("\\.rds$", "_linear.png", outrds), width = 8, height = 8,
    units = "in", res = 300)

print(ggplot(combcov$jcovscaled, aes(x = scaledcoverage, y = uniqreads)) + 
        geom_point(size = 0.5) + facet_wrap(~ method) + theme_bw() + 
        geom_abline(slope = 1, intercept = 0) + xlab("Scaled predicted coverage") + 
        ylab("Observed junction coverage"))

dev.off()

png(gsub("\\.rds$", "_log.png", outrds), width = 8, height = 8,
    units = "in", res = 300)

print(ggplot(combcov$jcovscaled, aes(x = scaledcoverage + 1, y = uniqreads + 1)) + 
        geom_point(size = 0.5) + facet_wrap(~ method) + theme_bw() + 
        geom_abline(slope = 1, intercept = 0) + xlab("Scaled predicted coverage") + 
        ylab("Observed junction coverage") + scale_y_log10() + scale_x_log10())

dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
