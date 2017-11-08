bf <- "/home/Shared/data/seq/hussain_bath_nanopore_rnaseq/NSK007/minimap2txome/SS2_wt_1/SS2_wt_1_minimap_txome_s.bam"

suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(dplyr))

reads <- readGAlignments(bf, use.names = TRUE)

alltx <- unique(as.character(seqnames(reads)))

uniqreads <- structure(rep(0, length(alltx)), names = alltx)
mmreads <- structure(rep(0, length(alltx)), names = alltx)

df <- data.frame(read = names(reads),
                 tx = as.character(seqnames(reads)),
                 stringsAsFactors = FALSE) %>%
  dplyr::distinct()

for (i in unique(df$read)) {
  tmp <- subset(df, read == i)
  if (nrow(tmp) == 1) {
    uniqreads[tmp$tx] <- uniqreads[tmp$tx] + 1
  } else {
    mmreads[tmp$tx] <- mmreads[tmp$tx] + 1
  }
}

allreads <- as.data.frame(uniqreads) %>% tibble::rownames_to_column("transcript") %>%
  dplyr::full_join(as.data.frame(mmreads) %>% tibble::rownames_to_column("transcript"))
allreads$transcript <- gsub("\\.[0-9]+", "", allreads$transcript)
