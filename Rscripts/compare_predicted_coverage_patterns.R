args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(predcov1)
print(predcov2)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
})

predcov1 <- readRDS(predcov1)
predcov2 <- readRDS(predcov2)

transcripts <- intersect(names(predcov1), names(predcov2))

## Calculate correlation between coverage patterns
cov_corrs <- sapply(transcripts, function(tx) {
  c1 <- predcov1[[tx]]$pred.cov
  c2 <- predcov2[[tx]]$pred.cov
  if (!is.null(c1) && !(is.null(c2)))
    cor(c1, c2)
  else
    NA
})

## Fraction NAs
mean(is.na(cov_corrs))

## Plot
df <- data.frame(transcript = names(cov_corrs),
                 cov_corrs = cov_corrs,
                 stringsAsFactors = FALSE)

pdf(gsub("rds$", "pdf", outrds))
print(ggplot(df, aes(x = cov_corrs)) + geom_histogram(bins = 100, fill = "lightblue") + 
        theme_bw() + xlab("Correlation of coverage patterns between samples") + 
        ylab("Frequency"))
dev.off()

saveRDS(df, file = outrds)
date()
sessionInfo()
