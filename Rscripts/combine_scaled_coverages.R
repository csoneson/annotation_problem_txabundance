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
print(genecharacteristics) ## Gene-level characteristics
print(exoncountstxt) ## Exon counts for each gene
print(introncountstxt) ## Intron counts for each gene
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
})

## Combine scaled junction coverages for all methods. Also summarize transcript
## and gene abundances, as well as characteristics of these features.

## Read junction coverages from STAR. Both in a stranded and unstranded mode.
jcov <- read.delim(junctioncovSTAR, header = FALSE, as.is = TRUE)
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

## Count the number of junctions per gene
jcovscaled <- jcovscaled %>% 
  dplyr::group_by(gene) %>% 
  dplyr::mutate(nbr_junctions_in_gene = length(start))

jcovscaled <- jcovscaled %>%
  dplyr::group_by(seqnames, start, end, gene) %>%
  dplyr::mutate(transcript = paste(unique(strsplit(paste(unique(transcript), collapse = ","), 
                                                   ",")[[1]]), collapse = ",")) %>% ## to make sure we have the same transcript combination for a given junction across all methods
  dplyr::ungroup() %>% 
  dplyr::mutate(pred.cov = replace(pred.cov, is.na(pred.cov), 0)) %>%
  dplyr::left_join(jcov, by = c("seqnames", "start", "end", "strand")) %>%
  dplyr::mutate(uniqreads = replace(uniqreads, is.na(uniqreads), 0),
                mmreads = replace(mmreads, is.na(mmreads), 0)) %>%
  dplyr::mutate(fracunique = uniqreads/(uniqreads + mmreads)) %>% 
  dplyr::mutate(fracunique = replace(fracunique, is.na(fracunique), 1)) %>%
  dplyr::ungroup()

## Define a junction ID for each junction
j0 <- jcovscaled %>% dplyr::select(seqnames, start, end, gene) %>%
  dplyr::distinct() %>% dplyr::group_by(gene) %>% dplyr::arrange(start) %>%
  dplyr::mutate(junctionid = paste0("J", seq_len(length(start)))) %>%
  dplyr::ungroup()

jcovscaled <- jcovscaled %>% dplyr::left_join(j0) %>%
  dplyr::select(junctionid, everything(), transcript)

## Read and combine transcript quantifications
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

## Summarize abundances on gene level
allquants_gene <- allquants %>% dplyr::group_by(gene, method) %>%
  dplyr::mutate(TPMrel = TPM/sum(TPM)) %>%
  dplyr::mutate(TPMrel = replace(TPMrel, is.na(TPMrel), 0)) %>%
  dplyr::summarize(count = sum(count),
                   TPM = sum(TPM),
                   nbr_expressed_transcripts = sum(TPM > 0),
                   nbr_expressed_transcripts_5p = sum(TPMrel > 0.05)) %>%
  dplyr::ungroup()

## Add gene characteristics
genechars <- readRDS(genecharacteristics)
allquants_gene <- dplyr::left_join(allquants_gene, genechars, 
                                   by = c("gene" = "gene_id"))

## Add exon and intron counts
exoncounts <- read.delim(exoncountstxt, skip = 1, header = TRUE, as.is = TRUE) %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
  setNames(c("gene", "exoncount"))
introncounts <- read.delim(introncountstxt, skip = 1, header = TRUE, as.is = TRUE) %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
  setNames(c("gene", "introncount"))
allquants_gene <- dplyr::left_join(allquants_gene, 
                                   dplyr::full_join(exoncounts, introncounts)) %>% 
  dplyr::mutate(exoncount = replace(exoncount, is.na(exoncount), 0)) %>% 
  dplyr::mutate(introncount = replace(introncount, is.na(introncount), 0)) %>%
  dplyr::mutate(intron_exon_ratio = introncount/exoncount) %>%
  dplyr::mutate(intron_exon_ratio = replace(intron_exon_ratio, exoncount==0 & introncount==0, 0))

## Add total unique and multimapping junction reads per gene
totjunctionreads <- jcovscaled %>% dplyr::filter(method == "Salmon") %>%
  dplyr::select(gene, uniqreads, mmreads) %>%
  dplyr::group_by(gene) %>% dplyr::summarize(uniqjuncreads = sum(uniqreads),
                                             mmjuncreads = sum(mmreads))
allquants_gene <- dplyr::left_join(allquants_gene, totjunctionreads, by = "gene") %>%
  dplyr::mutate(uniqjuncreads = replace(uniqjuncreads, is.na(uniqjuncreads), 0)) %>%
  dplyr::mutate(mmjuncreads = replace(mmjuncreads, is.na(mmjuncreads), 0)) %>%
  dplyr::mutate(uniqjuncfraction = uniqjuncreads/(uniqjuncreads + mmjuncreads)) %>%
  dplyr::mutate(uniqjuncfraction = replace(uniqjuncfraction, is.na(uniqjuncfraction), 1))

## Add number of junctions per gene (from the junction summary). If a gene is
## not in the junction summary, it means that it doesn't have any annotated
## junctions.
allquants_gene <- allquants_gene %>% 
  dplyr::left_join(jcovscaled %>% dplyr::select(gene, method, nbr_junctions_in_gene) %>%
                     dplyr::distinct(),
                   by = c("gene", "method")) %>%
  dplyr::mutate(nbr_junctions_in_gene = replace(nbr_junctions_in_gene, is.na(nbr_junctions_in_gene), 0))

saveRDS(list(junctions = jcovscaled, transcripts = allquants, 
             genes = allquants_gene), file = outrds)

sessionInfo()
date()
