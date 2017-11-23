args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(outrds)

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

gtf <- rtracklayer::import(gtf, format = "gtf")

## Extract only unique 3' UTRs by gene
utrs <- as.data.frame(gtf) %>% dplyr::select(seqnames, start, end, strand, gene_id, type) %>%
  dplyr::filter(type == "three_prime_utr") %>% dplyr::distinct()

## Summarize by gene
utrsum <- utrs %>% dplyr::group_by(gene_id) %>% dplyr::mutate(width = end - start + 1) %>%
  dplyr::summarize(lengthdiff = max(width) - min(width), 
                   nUTR = length(type))

pdf(gsub("rds$", "pdf", outrds))

## Number of distinct 3' UTRs per gene
ggplot(utrsum, aes(x = nUTR)) + geom_histogram(bins = 100) + theme_bw() + 
  xlab("Number of distinct 3' UTRs per gene")
ggplot(utrsum, aes(x = nUTR)) + geom_histogram(bins = 100) + theme_bw() + scale_y_sqrt() + 
  xlab("Number of distinct 3' UTRs per gene")
ggplot(utrsum %>% dplyr::mutate(nUTRgr = Hmisc::cut2(nUTR, g = 10)), aes(x = nUTRgr)) +
  geom_bar() + theme_bw() + 
  xlab("Number of distinct 3' UTRs per gene")

## Difference between longest and shortest 3' UTR
ggplot(utrsum, aes(x = lengthdiff)) + geom_histogram(bins = 100) + theme_bw() + 
  scale_y_sqrt() + xlab("Length difference between longest and shortest 3' UTR")
ggplot(utrsum %>% dplyr::mutate(lengthdiffgr = Hmisc::cut2(lengthdiff, g = 10)),
       aes(x = lengthdiffgr)) + geom_bar() + theme_bw() + 
  xlab("Length difference between longest and shortest 3' UTR")

dev.off()

sessionInfo()
date()
