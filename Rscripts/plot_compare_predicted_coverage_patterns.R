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

## Check whether any of the coverage patterns is uniform (if the prediction failed)
cov_types <- sapply(transcripts, function(tx) {
  c1 <- predcov1[[tx]]$note
  c2 <- predcov2[[tx]]$note
  if (c1 == "covOK" && c2 == "covOK")
    "noneUniform"
  else if (all(c(c1, c2) != "covOK"))
    "bothUniform"
  else 
    "oneUniform"
})

## Fraction NA correlations
mean(is.na(cov_corrs))
table(cov_types)

## Plot
stopifnot(all(names(cov_corrs) == names(cov_types)))
df <- data.frame(transcript = names(cov_corrs),
                 cov_corrs = cov_corrs,
                 cov_types = cov_types,
                 stringsAsFactors = FALSE)

pdf(gsub("rds$", "pdf", outrds))
print(ggplot(df, aes(x = cov_corrs)) + geom_histogram(bins = 100, fill = "lightblue") + 
        theme_bw() + xlab("Correlation of coverage patterns between samples") + 
        ylab("Frequency"))
print(ggplot(df, aes(x = cov_corrs)) + geom_histogram(bins = 100, fill = "lightblue") + 
        theme_bw() + xlab("Correlation of coverage patterns between samples") + 
        ylab("Frequency") + facet_wrap(~cov_types, ncol = 1))
dev.off()

saveRDS(df, file = outrds)
date()
sessionInfo()
