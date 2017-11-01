args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(junctioncov)
print(quantsf)
print(outrds)

suppressPackageStartupMessages(library(dplyr))
source("Rscripts/plot_tracks.R")

## Create gene models for Gviz visualization
options(ucscChromosomeNames = FALSE)
genemodels_exon <- create_genemodels(gtf, seltype = "exon")
genemodels_cds <- create_genemodels(gtf, seltype = "CDS")

jcov <- read.delim(junctioncov, 
                   header = FALSE, as.is = TRUE)
colnames(jcov) <- c("seqnames", "start", "end", "strand", "motif", "annot", 
                    "uniqreads", "mmreads", "maxoverhang")
jcov <- jcov %>% dplyr::mutate(strand = replace(strand, strand == 1, "+")) %>%
  dplyr::mutate(strand = replace(strand, strand == 2, "-")) %>%
  dplyr::select(-motif, -annot, -maxoverhang)

## Read Salmon quantifications
quantsf <- read.delim(quantsf, header = TRUE, as.is = TRUE)
quantsf$Name <- gsub("\\.[0-9]+", "", quantsf$Name)

saveRDS(list(genemodels_exon = genemodels_exon, genemodels_cds = genemodels_cds,
             jcov = jcov, quantsf = quantsf), file = outrds)

sessionInfo()
date()
