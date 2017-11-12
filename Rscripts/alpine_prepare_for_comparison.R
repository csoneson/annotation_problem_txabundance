args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(junctioncov) ## Junction reads from STAR
print(quantsf) ## Salmon quantifications
print(quantsfbwa) ## Salmon quantifications in alignment mode (after BWA)
print(heratsv) ## Hera quantifications
print(abundancetsv) ## kallisto quantifications
print(isoformsresults) ## RSEM quantifications
print(stringtiegtf) ## StringTie quantifications
print(quantsfnanopore) ## If present, nanopore results
print(outrds)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rtracklayer))
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

## Read effective lengths from Salmon (for StringTie later)
effl <- read.delim(quantsf, header = TRUE, as.is = TRUE) %>%
  dplyr::mutate(Name = gsub("\\.[0-9]+", "", Name)) %>%
  dplyr::rename(transcript = Name)

## Read nanopore quantifications
if (quantsfnanopore != "") {
  quantsfnanopore <- read.delim(quantsfnanopore, header = TRUE, as.is = TRUE) %>%
    dplyr::mutate(Name = gsub("\\.[0-9]+", "", Name)) %>%
    dplyr::select(Name, TPM, NumReads) %>%
    dplyr::rename(transcript = Name, count = NumReads) %>%
    dplyr::mutate(method = "SalmonNanopore")
} else {
  quantsfnanopore <- NULL
}

## Read Salmon quantifications
quantsf <- read.delim(quantsf, header = TRUE, as.is = TRUE) %>%
  dplyr::mutate(Name = gsub("\\.[0-9]+", "", Name)) %>%
  dplyr::select(Name, TPM, NumReads) %>%
  dplyr::rename(transcript = Name, count = NumReads) %>%
  dplyr::mutate(method = "Salmon")

## Read Salmon+BWA quantifications
quantsfbwa <- read.delim(quantsfbwa, header = TRUE, as.is = TRUE) %>%
  dplyr::mutate(Name = gsub("\\.[0-9]+", "", Name)) %>%
  dplyr::select(Name, TPM, NumReads) %>%
  dplyr::rename(transcript = Name, count = NumReads) %>%
  dplyr::mutate(method = "SalmonBWA")

## Read kallisto quantifications
abundancetsv <- read.delim(abundancetsv, header = TRUE, as.is = TRUE) %>%
  dplyr::mutate(target_id = gsub("\\.[0-9]+", "", target_id)) %>%
  dplyr::select(target_id, tpm, est_counts) %>%
  dplyr::rename(transcript = target_id, count = est_counts, TPM = tpm) %>%
  dplyr::mutate(method = "kallisto")

## Read hera quantifications
heratsv <- read.delim(heratsv, header = TRUE, as.is = TRUE) %>%
  dplyr::mutate(target_id = sapply(strsplit(X.target_id, ":"), .subset, 1)) %>%
  dplyr::select(target_id, tpm, est_counts) %>%
  dplyr::rename(transcript = target_id, count = est_counts, TPM = tpm) %>%
  dplyr::mutate(method = "Hera")

## Read RSEM quantifications
isoformsresults <- read.delim(isoformsresults, header = TRUE, as.is = TRUE) %>%
  dplyr::mutate(transcript_id = gsub("\\.[0-9]+", "", transcript_id)) %>%
  dplyr::select(transcript_id, TPM, expected_count) %>%
  dplyr::rename(transcript = transcript_id, count = expected_count) %>%
  dplyr::mutate(method = "RSEM")

## Read StringTie quantifications
## Use Salmon effective lengths to derive an "expected count" for StringTie
stringtiegtf <- import(stringtiegtf) 
stringtiegtf <- subset(stringtiegtf, type == "transcript")
stringtiegtf <- as.data.frame(mcols(stringtiegtf)[, c("transcript_id", "TPM")]) %>%
  dplyr::rename(transcript = transcript_id) %>%
  dplyr::mutate(TPM = as.numeric(as.character(TPM))) %>%
  dplyr::left_join(effl %>% dplyr::select(transcript, EffectiveLength)) %>%
  dplyr::mutate(count = TPM * EffectiveLength/sum(TPM * EffectiveLength, na.rm = TRUE) * sum(effl$NumReads)) %>%
  dplyr::select(transcript, TPM, count) %>%
  dplyr::mutate(method = "StringTie")

## Make a list with one element for each method, with columns transcript, TPM, count
quants <- rbind(quantsf, quantsfbwa, abundancetsv, heratsv, isoformsresults, 
                stringtiegtf, quantsfnanopore)

saveRDS(list(genemodels_exon = genemodels_exon, genemodels_cds = genemodels_cds,
             jcov = jcov, quants = quants), file = outrds)

sessionInfo()
date()
