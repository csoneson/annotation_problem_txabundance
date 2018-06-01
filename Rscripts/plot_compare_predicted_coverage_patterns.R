################################################################################
##                                                                            ##
## Compare coverage patterns between samples                                  ##
##                                                                            ##
## Inputs:                                                                    ##
## * predcov1: object with predicted transcript coverage patterns for first   ##
##             sample                                                         ##
## * predcov2: object with predicted transcript coverage patterns for second  ##
##             sample                                                         ##
## * samplename1: the sample name for the predcov1 file                       ##
## * samplename2: the sample name for the predcov2 file                       ##
## * outrds: output file                                                      ##
##                                                                            ##
## Outputs:                                                                   ##
## * png figures with the distribution of correlation coefficients for        ##
##   transcript coverage patterns in the two samples, across all transcripts  ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(predcov1)
print(predcov2)
print(samplename1)
print(samplename2)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
  library(cowplot)
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

## Get transcripts with highest and lowest correlation
txhigh <- df$transcript[which.max(df$cov_corrs)]
txlow <- df$transcript[which.min(df$cov_corrs)]
df1 <- rbind(data.frame(x = seq_len(length(predcov1[[txhigh]]$pred.cov)),
                        Coverage = predcov1[[txhigh]]$pred.cov/max(predcov1[[txhigh]]$pred.cov),
                        sample = samplename1,
                        dtype = "High correlation",
                        stringsAsFactors = FALSE),
             data.frame(x = seq_len(length(predcov2[[txhigh]]$pred.cov)),
                        Coverage = predcov2[[txhigh]]$pred.cov/max(predcov2[[txhigh]]$pred.cov),
                        sample = samplename2, 
                        dtype = "High correlation",
                        stringsAsFactors = FALSE),
             data.frame(x = seq_len(length(predcov1[[txlow]]$pred.cov)),
                        Coverage = predcov1[[txlow]]$pred.cov/max(predcov1[[txlow]]$pred.cov),
                        sample = samplename1,
                        dtype = "Low correlation",
                        stringsAsFactors = FALSE),
             data.frame(x = seq_len(length(predcov2[[txlow]]$pred.cov)),
                        Coverage = predcov2[[txlow]]$pred.cov/max(predcov2[[txlow]]$pred.cov),
                        sample = samplename2, 
                        dtype = "Low correlation",
                        stringsAsFactors = FALSE))

png(gsub("\\.rds$", "_all.png", outrds), height = 6, width = 9, unit = "in", res = 300)
print(plot_grid(
  ggplot(df, aes(x = cov_corrs)) + geom_histogram(bins = 100, fill = "#7BAFDE") + 
    theme_bw() + xlab("Correlation of coverage patterns between samples") + 
    ylab("Frequency"),
  ggplot(df1, aes(x = x, y = Coverage, color = sample)) + 
    geom_line(size = 1.25) + ylab("Predicted coverage") + 
    theme_bw() + xlab("Position in transcript") + 
    facet_wrap(~ dtype, ncol = 1, scales = "free") + 
    scale_color_manual(values = c("#882E72", "#90C987"), name = "") + 
    theme(legend.position = "bottom"),
  nrow = 1, rel_widths = c(2, 1), labels = c("A", "B")
))
dev.off()

png(gsub("\\.rds$", "_strat.png", outrds), height = 6, width = 6, unit = "in", res = 300)
print(ggplot(df, aes(x = cov_corrs)) + geom_histogram(bins = 100, fill = "lightblue") + 
        theme_bw() + xlab("Correlation of coverage patterns between samples") + 
        ylab("Frequency") + facet_wrap(~cov_types, ncol = 1))
dev.off()

saveRDS(df, file = outrds)
date()
sessionInfo()
