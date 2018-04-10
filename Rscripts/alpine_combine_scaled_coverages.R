args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

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
suppressPackageStartupMessages(library(ggplot2))

jcov <- read.delim(junctioncovSTAR, 
                   header = FALSE, as.is = TRUE)
colnames(jcov) <- c("seqnames", "start", "end", "strand", "motif", "annot", 
                    "uniqreads", "mmreads", "maxoverhang")
jcov <- jcov %>% dplyr::mutate(strand = replace(strand, strand == 1, "+")) %>%
  dplyr::mutate(strand = replace(strand, strand == 2, "-")) %>%
  dplyr::select(seqnames, start, end, strand, uniqreads, mmreads)

jcovnostrand <- jcov %>% group_by(seqnames, start, end) %>%
  dplyr::summarize(uniqreads = sum(uniqreads), 
                   mmreads = sum(mmreads)) %>%
  dplyr::mutate(strand = "*") %>% dplyr::ungroup() %>%
  dplyr::select(seqnames, start, end, strand, uniqreads, mmreads) %>%
  as.data.frame()
jcov <- rbind(jcov, jcovnostrand)

## Read and merge junction coverages predicted by each method
jcovscaled <- do.call(rbind, list(readRDS(junctioncovSalmon)$allcovs,
                                  readRDS(junctioncovSalmonBWA)$allcovs,
                                  readRDS(junctioncovhera)$allcovs,
                                  readRDS(junctioncovkallisto)$allcovs,
                                  readRDS(junctioncovRSEM)$allcovs,
                                  readRDS(junctioncovStringTie)$allcovs))
if (junctioncovNanopore != "") {
  jcovscaled <- rbind(jcovscaled, readRDS(junctioncovNanopore)$allcovs)
}
if (junctioncovSalmonCDS != "") {
  jcovscaled <- rbind(jcovscaled, readRDS(junctioncovSalmonCDS)$allcovs)
}

jcovscaled <- jcovscaled %>%
  dplyr::group_by(seqnames, start, end, gene) %>%
  dplyr::mutate(transcript = paste(unique(strsplit(paste(unique(transcript), collapse = ","), 
                                                   ",")[[1]]), collapse = ",")) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(pred.cov = replace(pred.cov, is.na(pred.cov), 0)) %>%
  dplyr::left_join(jcov, by = c("seqnames", "start", "end", "strand")) %>%
  dplyr::mutate(uniqreads = replace(uniqreads, is.na(uniqreads), 0),
                mmreads = replace(mmreads, is.na(mmreads), 0)) %>%
  dplyr::group_by(gene, method) %>% 
  dplyr::mutate(scaledcoverage = pred.cov/sum(pred.cov, na.rm = TRUE) * 
                  sum(uniqreads, na.rm = TRUE)) %>%
  dplyr::mutate(scaledcoverage = replace(scaledcoverage, is.na(scaledcoverage), 0)) %>% 
  dplyr::mutate(score = round(sum(abs(uniqreads - scaledcoverage), na.rm = TRUE)/
                                sum(uniqreads, na.rm = TRUE), 2)) %>% 
  dplyr::mutate(methodscore = paste0(method, " (", score, ")")) %>%
  dplyr::ungroup()

j0 <- jcovscaled %>% dplyr::select(seqnames, start, end, gene) %>%
  dplyr::distinct() %>% dplyr::group_by(gene) %>% dplyr::arrange(start) %>%
  dplyr::mutate(junctionid = paste0("J", seq_len(length(start)))) %>%
  dplyr::ungroup()

jcovscaled <- jcovscaled %>% dplyr::left_join(j0) %>%
  dplyr::select(junctionid, everything(), transcript)

allquants <- do.call(rbind, list(readRDS(junctioncovSalmon)$quants,
                                 readRDS(junctioncovSalmonBWA)$quants,
                                 readRDS(junctioncovhera)$quants,
                                 readRDS(junctioncovkallisto)$quants,
                                 readRDS(junctioncovRSEM)$quants,
                                 readRDS(junctioncovStringTie)$quants))

if (junctioncovNanopore != "") {
  allquants <- rbind(allquants, readRDS(junctioncovNanopore)$quants)
}
if (junctioncovSalmonCDS != "") {
  allquants <- rbind(allquants, readRDS(junctioncovSalmonCDS)$quants)
}

saveRDS(list(jcov = jcov, jcovscaled = jcovscaled, allquants = allquants), file = outrds)

sessionInfo()
date()
