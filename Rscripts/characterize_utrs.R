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
utrs <- as.data.frame(gtf) %>% dplyr::filter(type == "three_prime_utr") %>%
  dplyr::mutate(width = end - start + 1) %>% dplyr::group_by(transcript_id) %>%
  dplyr::summarize(start = min(start), end = max(start), width = sum(width), 
                   seqnames = seqnames[1], gene_id = gene_id[1]) %>%
  dplyr::select(-transcript_id) %>% 
  dplyr::group_by(gene_id) %>% dplyr::summarize(lengthdiff = max(width) - min(width),
                                                nUTR = length(width))

pdf(gsub("rds$", "pdf", outrds))

## Number of distinct 3' UTRs per gene
ggplot(utrs, aes(x = nUTR)) + geom_histogram(bins = 100) + theme_bw() + 
  xlab("Number of distinct 3' UTRs per gene")
ggplot(utrs, aes(x = nUTR)) + geom_histogram(bins = 100) + theme_bw() + scale_y_sqrt() + 
  xlab("Number of distinct 3' UTRs per gene")
ggplot(utrs %>% dplyr::mutate(nUTRgr = Hmisc::cut2(nUTR, g = 10)), aes(x = nUTRgr)) +
  geom_bar() + theme_bw() + 
  xlab("Number of distinct 3' UTRs per gene")

## Difference between longest and shortest 3' UTR
ggplot(utrs, aes(x = lengthdiff)) + geom_histogram(bins = 100) + theme_bw() + 
  scale_y_sqrt() + xlab("Length difference between longest and shortest 3' UTR")
ggplot(utrs %>% dplyr::mutate(lengthdiffgr = Hmisc::cut2(lengthdiff, g = 10)),
       aes(x = lengthdiffgr)) + geom_bar() + theme_bw() + 
  xlab("Length difference between longest and shortest 3' UTR")

dev.off()

saveRDS(utrs, file = outrds)
sessionInfo()
date()
