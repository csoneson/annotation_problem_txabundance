args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(geneids)
print(bwfile)
print(gtffile)
print(outrds)

geneids <- strsplit(geneids, ",")[[1]]

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Gviz))
source("Rscripts/plot_tracks.R")

options(ucscChromosomeNames = FALSE)

genemodels_exon <- create_genemodels(gtffile, seltype = "exon")
genemodels_cds <- create_genemodels(gtffile, seltype = "CDS")

pdf(gsub("rds$", "pdf", outrds), width = 12, height = 8)
for (gid in geneids) {
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

sessionInfo()
date()
