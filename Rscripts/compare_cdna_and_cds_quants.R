args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(cdnaquant)
print(cdsquant)
print(gtffile)
print(tx2gene)
print(bwfile)
print(outrds)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Gviz))
source("Rscripts/plot_tracks.R")
options(ucscChromosomeNames = FALSE)

## Read cDNA quantifications
cdna <- read.delim(cdnaquant, header = TRUE, as.is = TRUE)

## Read CDS quantifications
cds <- read.delim(cdsquant, header = TRUE, as.is = TRUE)

## Read tx2gene
tx2gene <- readRDS(tx2gene)

## Add gene assignment
cdna$gene <- tx2gene$gene[match(cdna$Name, tx2gene$tx)]
cds$gene <- tx2gene$gene[match(cds$Name, tx2gene$tx)]

## Calculate relative abundances
cdna <- cdna %>% dplyr::group_by(gene) %>% 
  dplyr::mutate(RelCount = NumReads/sum(NumReads))
cds <- cds %>% dplyr::group_by(gene) %>% 
  dplyr::mutate(RelCount = NumReads/sum(NumReads))

## Find shared transcripts
cdna <- cdna %>% dplyr::select(Name, gene, NumReads, RelCount) %>%
  dplyr::rename(NumReadscDNA = NumReads, RelCountcDNA = RelCount)
cds <- cds %>% dplyr::select(Name, gene, NumReads, RelCount) %>%
  dplyr::rename(NumReadsCDS = NumReads, RelCountCDS = RelCount)

abundances <- dplyr::full_join(cdna, cds) %>% 
  dplyr::mutate(MaxNumReads = pmax(NumReadscDNA, NumReadsCDS)) %>%
  dplyr::group_by(gene) %>% 
  dplyr::mutate(Ntx = length(gene)) %>%
  dplyr::ungroup()

tophits <- abundances %>% dplyr::filter(MaxNumReads > 2000) %>% 
  dplyr::filter(Ntx > 1) %>% 
  dplyr::filter(abs(RelCountcDNA - RelCountCDS) > 0.9)

write.table(tophits, file = gsub("rds$", "txt", outrds), row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = "\t")

## Plot coverage of top hits
(geneids <- unique(gsub("\\.[0-9]+$", "", tophits$gene)))

genemodels_exon <- create_genemodels(gtffile, seltype = "exon")
genemodels_cds <- create_genemodels(gtffile, seltype = "CDS")

pdf(gsub("rds$", "pdf", outrds), width = 12, height = 8)
for (gid in geneids) {
  message(gid)
  tryCatch({
    plot_tracks(mygene = gid, genemodels = genemodels_exon, 
                genemodels2 = genemodels_cds, 
                gtf_file = NULL, rnaseq_datafiles = structure(bwfile, names = "s1"), 
                rnaseq_condition = structure("g1", names = "s1"), show_chr = NULL, 
                min_coord = NULL, max_coord = NULL, 
                pdf_filename = NULL, pdf_width = 7, pdf_height = 7)
  }, error = function(e) message(e))
}
dev.off()

saveRDS(list(abundances = abundances, tophits = tophits), file = outrds)

sessionInfo()
date()
