args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(outrds)

source("Rscripts/helper_plot_tracks.R")

## Create gene models for Gviz visualization
options(ucscChromosomeNames = FALSE)
genemodels_exon <- create_genemodels(gtf, seltype = "exon")
genemodels_cds <- create_genemodels(gtf, seltype = "CDS")

saveRDS(list(genemodels_exon = genemodels_exon, genemodels_cds = genemodels_cds), 
        file = outrds)

sessionInfo()
date()
