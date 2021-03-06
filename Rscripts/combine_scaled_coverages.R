################################################################################
##                                                                            ##
## Combine scaled junction coverages from different methods                   ##
##                                                                            ##
## Inputs:                                                                    ##
## * junctioncovSTAR: junction coverages from STAR (uniquely and multimapping ##
##                    reads)                                                  ##
## * junctioncovSalmon: junction coverages and abundances from Salmon         ##
## * junctioncovSalmonSTAR: junction coverages and abundances from SalmonSTAR ##
## * junctioncovSalmonKeepDup: junction coverages and abundances from         ##
##                             SalmonKeepDup                                  ##
## * junctioncovSalmonCDS: junction coverages and abundances from SalmonCDS   ##
## * junctioncovhera: junction coverages and abundances from hera             ##
## * junctioncovkallisto: junction coverages and abundances from kallisto     ##
## * junctioncovRSEM: junction coverages and abundances from RSEM             ##
## * junctioncovStringTie: junction coverages and abundances from StringTie   ##
## * junctioncovNanopore: junction coverages and abundances from Nanopore     ##
## * genecharacteristics: table with gene characteristics                     ##
## * exoncountstxt: text file with exon counts for each gene                  ##
## * introncountstxt: text file with intron counts for each gene              ##
## * outrds: output file                                                      ##
##                                                                            ##
## Outputs:                                                                   ##
## * A table with combined junction coverages, one with all transcript        ##
##   abundances and one with gene abundance estimates                         ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(junctioncovSTAR)
print(junctioncovSalmon)
print(junctioncovSalmonSTAR)
print(junctioncovSalmonKeepDup)
print(junctioncovSalmonCDS)
print(junctioncovhera)
print(junctioncovkallisto)
print(junctioncovRSEM)
print(junctioncovStringTie)
print(junctioncovNanopore)
print(genecharacteristics)
print(exoncountstxt)
print(introncountstxt)
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
                                  readRDS(junctioncovkallisto)$allcovs))
if (junctioncovSalmonSTAR != "") {
  jcovscaled <- rbind(jcovscaled, readRDS(junctioncovSalmonSTAR)$allcovs)
}
if (junctioncovhera != "") {
  jcovscaled <- rbind(jcovscaled, readRDS(junctioncovhera)$allcovs)
}
if (junctioncovRSEM != "") {
  jcovscaled <- rbind(jcovscaled, readRDS(junctioncovRSEM)$allcovs)
}
if (junctioncovStringTie != "") {
  jcovscaled <- rbind(jcovscaled, readRDS(junctioncovStringTie)$allcovs)
}
if (junctioncovNanopore != "") {
  jcovscaled <- rbind(jcovscaled, readRDS(junctioncovNanopore)$allcovs)
}
if (junctioncovSalmonCDS != "") {
  jcovscaled <- rbind(jcovscaled, readRDS(junctioncovSalmonCDS)$allcovs)
}
if (junctioncovSalmonKeepDup != "") {
  jcovscaled <- rbind(jcovscaled, readRDS(junctioncovSalmonKeepDup)$allcovs)
}

## Count the number of junctions per gene
jcovscaled <- jcovscaled %>% 
  dplyr::group_by(gene, method) %>% 
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
                                 readRDS(junctioncovkallisto)$quants))

if (junctioncovSalmonSTAR != "") {
  allquants <- rbind(allquants, readRDS(junctioncovSalmonSTAR)$quants)
}
if (junctioncovhera != "") {
  allquants <- rbind(allquants, readRDS(junctioncovhera)$quants)
}
if (junctioncovRSEM != "") {
  allquants <- rbind(allquants, readRDS(junctioncovRSEM)$quants)
}
if (junctioncovStringTie != "") {
  allquants <- rbind(allquants, readRDS(junctioncovStringTie)$quants)
}
if (junctioncovNanopore != "") {
  allquants <- rbind(allquants, readRDS(junctioncovNanopore)$quants)
}
if (junctioncovSalmonCDS != "") {
  allquants <- rbind(allquants, readRDS(junctioncovSalmonCDS)$quants)
}
if (junctioncovSalmonKeepDup != "") {
  allquants <- rbind(allquants, readRDS(junctioncovSalmonKeepDup)$quants)
}

## Summarize abundances on gene level
allquants_gene <- allquants %>% dplyr::group_by(gene, method) %>%
  dplyr::mutate(TPMrel = TPM/sum(TPM)) %>%
  dplyr::mutate(TPMrel = replace(TPMrel, is.na(TPMrel), 0)) %>%
  dplyr::summarize(count = sum(count),
                   TPM = sum(TPM),
                   nbr_expressed_transcripts = sum(TPM > 0),
                   nbr_expressed_transcripts_5p = sum(TPMrel > 0.05),
                   covOKfraction = length(which(covnote == "covOK"))/length(covnote)) %>%
  dplyr::ungroup()

## Add gene characteristics
if (genecharacteristics != "") {
  genechars <- readRDS(genecharacteristics)
  allquants_gene <- dplyr::left_join(allquants_gene, genechars, 
                                     by = c("gene" = "gene_id"))
}

## Add exon and intron counts
if (exoncountstxt != "" && introncountstxt != "") {
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
}

## Add total unique and multimapping junction reads per gene
if ("Salmon" %in% jcovscaled$method) {
  selmethod <- "Salmon"
} else {
  selmethod <- "SalmonPermuted"
} 
totjunctionreads <- jcovscaled %>% dplyr::filter(method == selmethod) %>%
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
