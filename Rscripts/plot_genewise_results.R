args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gene)  ## gene of interest, or file listing collection of genes (one per row)
print(bigwig)  ## bigwig file for visualization
print(bigwignanopore)  ## second bigwig file (nanopore)
print(genemodels)  ## gene models 
print(combcovrds)  ## combined junction coverages
print(ncores)  ## number of cores for parallel computations
print(outdir)  ## output directory
print(libid)  ## string that will be added to the beginning of all output files (library ID)
print(checkdir)  ## directory to write (empty) rds files (time stamps)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
})

source("Rscripts/helper_plot_tracks.R")

## Read gene models for Gviz plot (pregenerated from gtf to save time) and
## quantifications
genemodels <- readRDS(genemodels)
combcov <- readRDS(combcovrds)

## Determine which gene(s) to investigate
if (file.exists(gene)) {
  genes <- unlist(read.delim(gene, as.is = TRUE, header = FALSE))
} else {
  genes <- gene
}

## Investigate each gene
mclapply(genes, function(currgene) {
  jl <- combcov$junctions %>% dplyr::filter(gene == currgene) %>%
    dplyr::mutate(junctionid2 = junctionid) %>%
    dplyr::mutate(junctionid2 = replace(junctionid2, 
                                        abs(scaled.cov - uniqreads) < mean(uniqreads), 
                                        ""))

  pdf(paste0(outdir, "/plots/", libid, currgene, ".pdf"), width = 12, height = 10)
  bwfiles <- structure(bigwig, names = "Illumina")
  bwcond <- structure("g1", names = "Illumina")
  if (bigwignanopore != "") {
    bwfiles <- c(bwfiles, structure(bigwignanopore, names = "Nanopore"))
    bwcond <- c(bwcond, structure("g2", names = "Nanopore"))
  }
  tryCatch({
    plot_tracks(mygene = currgene, genemodels = genemodels$genemodels_exon, 
                genemodels2 = genemodels$genemodels_cds, 
                gtf_file = NULL, rnaseq_datafiles = bwfiles, 
                rnaseq_condition = bwcond, show_chr = NULL, 
                min_coord = NULL, max_coord = NULL, 
                pdf_filename = NULL, pdf_width = 7, pdf_height = 7)
  }, error = function(e) message(e))
  
  print(ggplot(jl, aes(x = scaled.cov, y = uniqreads, label = junctionid)) + 
          geom_point(size = 4) + geom_label_repel() +
          facet_wrap(~ methodscore) + 
          geom_abline(intercept = 0, slope = 1) + 
          xlab("Scaled predicted coverage") + 
          ylab("Number of uniquely mapped reads") + 
          theme_bw())
  
  print(ggplot(jl, aes(x = scaled.cov, y = uniqreads, label = junctionid,
                       color = fracunique > 0.75)) + 
          geom_point(size = 4) + geom_label_repel() +
          facet_wrap(~ methodscore) + 
          geom_abline(intercept = 0, slope = 1) + 
          xlab("Scaled predicted coverage") + 
          ylab("Number of uniquely mapped reads") + 
          theme_bw())
  
  print(ggplot(jl, aes(x = scaled.cov, y = uniqreads, label = junctionid2)) + 
          geom_point(size = 4) + geom_label_repel() +
          facet_wrap(~ methodscore) + 
          geom_abline(intercept = 0, slope = 1) + 
          xlab("Scaled predicted coverage") + 
          ylab("Number of uniquely mapped reads") + 
          theme_bw())
  dev.off()
  
  write.table(jl %>% dplyr::select(-score, -pred.cov, -method, -junctionid2) %>%
                dplyr::mutate(scaled.cov = round(scaled.cov, 2)) %>% 
                tidyr::spread(methodscore, scaled.cov) %>%
                dplyr::arrange(start),
              file = paste0(outdir, "/jcov/", libid, currgene, "_jscaledcov.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  write.table(combcov$transcripts %>% dplyr::filter(gene == currgene) %>% 
                dplyr::select(transcript, method, TPM) %>%
                dplyr::mutate(TPM = round(TPM, 2)) %>% 
                tidyr::spread(method, TPM),
              file = paste0(outdir, "/tpm/", libid, currgene, "_tpm.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  write.table(combcov$transcripts %>% dplyr::filter(gene == currgene) %>% 
                dplyr::select(transcript, method, count) %>%
                dplyr::mutate(count = round(count, 2)) %>% 
                tidyr::spread(method, count),
              file = paste0(outdir, "/count/", libid, currgene, "_count.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  saveRDS(jl, paste0(checkdir, "/", currgene, ".rds"))
}, mc.preschedule = FALSE, mc.cores = ncores)

saveRDS(NULL, paste0(checkdir, "/", basename(gene), ".rds"))
sessionInfo()
date()

