args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(junctioncovSTAR) ## Junction reads from STAR
print(junctioncovSalmon) ## Salmon quantifications
print(junctioncovSalmonBWA) ## Salmon quantifications in alignment mode (after BWA)
print(junctioncovSalmonCDS) ## Salmon quantifications based on CDSs only
print(junctioncovhera) ## Hera quantifications
print(junctioncovkallisto) ## kallisto quantifications
print(junctioncovRSEM) ## RSEM quantifications
print(junctioncovStringTie) ## StringTie quantifications
print(junctioncovNanopore) ## If present, nanopore results
print(outrds)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rtracklayer))
source("Rscripts/plot_tracks.R")

## Create gene models for Gviz visualization
options(ucscChromosomeNames = FALSE)
genemodels_exon <- create_genemodels(gtf, seltype = "exon")
genemodels_cds <- create_genemodels(gtf, seltype = "CDS")

jcov <- read.delim(junctioncovSTAR, 
                   header = FALSE, as.is = TRUE)
colnames(jcov) <- c("seqnames", "start", "end", "strand", "motif", "annot", 
                    "uniqreads", "mmreads", "maxoverhang")
jcov <- jcov %>% dplyr::mutate(strand = replace(strand, strand == 1, "+")) %>%
  dplyr::mutate(strand = replace(strand, strand == 2, "-")) %>%
  dplyr::select(-motif, -annot, -maxoverhang)

## Read and merge junction coverages predicted by each method
jcovscaled <- do.call(rbind, list(readRDS(junctioncovSalmon)$allcovs,
                                  readRDS(junctioncovSalmonBWA)$allcovs,
                                  readRDS(junctioncovSalmonCDS)$allcovs,
                                  readRDS(junctioncovhera)$allcovs,
                                  readRDS(junctioncovkallisto)$allcovs,
                                  readRDS(junctioncovRSEM)$allcovs,
                                  readRDS(junctioncovStringTie)$allcovs))
if (junctioncovNanopore != "") {
  jcovscaled <- rbind(jcovscaled, readRDS(junctioncovNanopore)$allcovs)
}

allquants <- do.call(rbind, list(readRDS(junctioncovSalmon)$quants,
                                 readRDS(junctioncovSalmonBWA)$quants,
                                 readRDS(junctioncovSalmonCDS)$quants, 
                                 readRDS(junctioncovhera)$quants,
                                 readRDS(junctioncovkallisto)$quants,
                                 readRDS(junctioncovRSEM)$quants,
                                 readRDS(junctioncovStringTie)$quants))

if (junctioncovNanopore != "") {
  allquants <- rbind(allquants, readRDS(junctioncovNanopore)$quants)
}

saveRDS(list(genemodels_exon = genemodels_exon, genemodels_cds = genemodels_cds,
             jcov = jcov, jcovscaled = jcovscaled, allquants = allquants), file = outrds)

sessionInfo()
date()
