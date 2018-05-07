args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(covrds)
print(outtxt)  ## output file

suppressPackageStartupMessages({
  library(dplyr)
})

covr <- readRDS(covrds)

covprednote <- sapply(covr, function(w) w$note)

write.table(
  dplyr::full_join(
    data.frame(table(covprednote)),
    data.frame(covprednote = c("covError", "covNA", "covOK"),
               explanation = c("transcripts for which the coverage prediction gave an error, e.g. if the transcript is shorter than the fragment length", 
                               "transcripts for which the coverage prediction returned NULL, e.g. if there are no reads mapping in the region", 
                               "transcripts where the coverage could be predicted without issues")),
    by = "covprednote"
  ),
  file = outtxt, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

date()
sessionInfo()
