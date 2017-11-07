jbt <- readRDS("reference/junctions_by_transcript.rds")
bf <- "/home/Shared/data/seq/hussain_bath_nanopore_rnaseq/NSK007/minimap2genome/SS2_wt_1/SS2_wt_1_minimap_genome_s.bam"

suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(dplyr))

print(outrds)

x <- GenomicAlignments::readGAlignments(bf, use.names = TRUE, 
                                        param = Rsamtools::ScanBamParam(tag = "NH"))
junc <- GenomicAlignments::junctions(x)
M <- mclapply(junc, function(rj) { ## For each read
  m <- c()
  ## First crude search
  f <- GenomicRanges::findOverlaps(query = rj, subject = jbt, type = "any",
                                   ignore.strand = TRUE)
  ## More refined search
  for (sh in unique(subjectHits(f))) {
    if (length(unique(queryHits(f)[subjectHits(f) == sh])) == nLnode(f)) {
      fp <- GenomicRanges::findOverlaps(query = rj, subject = jbt[[sh]], 
                                        type = "equal", ignore.strand = TRUE,
                                        maxgap = 5)
      if (length(unique(queryHits(fp))) == nLnode(fp)) {
        m <- c(m, names(jbt)[sh])
      }
    }
  }
  m
}, mc.preschedule = FALSE, mc.cores = 10)

alltx <- unique(unlist(M))
uniqreads <- structure(rep(0, length(alltx)), names = alltx)
mmreads <- structure(rep(0, length(alltx)), names = alltx)
for (i in seq_along(M)) {
  if (!is.null(M[[i]])) {
    if (length(M[[i]]) == 1) {
      uniqreads[M[[i]]] <- uniqreads[M[[i]]] + 1
    }
    else {
      mmreads[M[[i]]] <- mmreads[M[[i]]] + 1
    }
  }
}

allreads <- as.data.frame(uniqreads) %>% tibble::rownames_to_column("transcript") %>%
  dplyr::full_join(as.data.frame(mmreads) %>% tibble::rownames_to_column("transcript"))

saveRDS(list(M = M, allreads = allreads), file = outrds)

sessionInfo()
date()
