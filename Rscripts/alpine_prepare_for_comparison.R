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
suppressPackageStartupMessages(library(ggplot2))
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
                                  readRDS(junctioncovSalmonCDS)$allcovs,
                                  readRDS(junctioncovhera)$allcovs,
                                  readRDS(junctioncovkallisto)$allcovs,
                                  readRDS(junctioncovRSEM)$allcovs,
                                  readRDS(junctioncovStringTie)$allcovs))
if (junctioncovNanopore != "") {
  jcovscaled <- rbind(jcovscaled, readRDS(junctioncovNanopore)$allcovs)
}

jcovscaled <- jcovscaled %>%
  dplyr::group_by(seqnames, start, end, gene) %>%
  dplyr::mutate(transcript = paste(unique(strsplit(paste(unique(transcript), collapse = ","), ",")[[1]]), collapse = ",")) %>% 
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
                                 readRDS(junctioncovSalmonCDS)$quants, 
                                 readRDS(junctioncovhera)$quants,
                                 readRDS(junctioncovkallisto)$quants,
                                 readRDS(junctioncovRSEM)$quants,
                                 readRDS(junctioncovStringTie)$quants))

if (junctioncovNanopore != "") {
  allquants <- rbind(allquants, readRDS(junctioncovNanopore)$quants)
}

pdf(gsub("rds$", "pdf", outrds))
## Score distribution
print(ggplot(jcovscaled %>% dplyr::select(gene, method, score) %>% dplyr::distinct(), 
             aes(x = score, color = method)) + geom_density() + theme_bw() + 
        ggtitle("Score distribution") + 
        scale_color_manual(values = c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                                      "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                                      "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                                      "#90C987")[seq_len(length(unique(jcovscaled$method)))]))
print(ggplot(jcovscaled %>% dplyr::select(gene, method, score) %>% dplyr::distinct(), 
             aes(x = method, y = score, color = method)) + geom_boxplot() + theme_bw() + 
        ggtitle("Score distribution") + 
        scale_color_manual(values = c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                                      "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                                      "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                                      "#90C987")[seq_len(length(unique(jcovscaled$method)))]))
## Rank distribution (1 is best)
print(ggplot(jcovscaled %>% dplyr::select(gene, method, score) %>% dplyr::distinct() %>% 
               dplyr::mutate(score = replace(score, is.na(score), 10)) %>% 
               dplyr::group_by(gene) %>% dplyr::mutate(rank = rank(score)) %>%
               dplyr::mutate(keep = !(var(rank) == 0)) %>% dplyr::filter(keep), 
             aes(x = rank)) + geom_bar() + facet_wrap(~method) + theme_bw())
dev.off()

saveRDS(list(genemodels_exon = genemodels_exon, genemodels_cds = genemodels_cds,
             jcov = jcov, jcovscaled = jcovscaled, allquants = allquants), file = outrds)

sessionInfo()
date()
