args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gene)  ## gene of interest, or file listing collection of genes (one per row)
print(bigwig)  ## bigwig file for visualization
print(genemodels)  ## gene models 
print(combcovrds)  ## combined junction coverages
print(ncores)  ## number of cores for parallel computations
print(outdir)  ## output directory
print(checkdir)  ## directory to write (empty) rds files (time stamps)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))

source("Rscripts/plot_tracks.R")

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
  jl <- combcov$jcovscaled %>% dplyr::filter(gene == currgene)

  pdf(paste0(outdir, "/plots/", currgene, ".pdf"), width = 12, height = 10)
  tryCatch({
    plot_tracks(mygene = currgene, genemodels = genemodels$genemodels_exon, 
                genemodels2 = genemodels$genemodels_cds, 
                gtf_file = NULL, rnaseq_datafiles = structure(bigwig, names = "s1"), 
                rnaseq_condition = structure("g1", names = "s1"), show_chr = NULL, 
                min_coord = NULL, max_coord = NULL, 
                pdf_filename = NULL, pdf_width = 7, pdf_height = 7)
  }, error = function(e) message(e))
  
  print(ggplot(jl, aes(x = scaledcoverage, y = uniqreads, label = junctionid)) + 
          geom_point(size = 4) + geom_label_repel() +
          facet_wrap(~ methodscore) + 
          geom_abline(intercept = 0, slope = 1) + 
          xlab("Scaled predicted coverage") + 
          ylab("Number of uniquely mapped reads") + 
          theme_bw())
  dev.off()
  
  write.table(jl %>% dplyr::select(-score, -pred.cov, -method) %>%
                dplyr::mutate(scaledcoverage = round(scaledcoverage, 2)) %>% 
                tidyr::spread(methodscore, scaledcoverage) %>%
                dplyr::arrange(start),
              file = paste0(outdir, "/jcov/", currgene, "_jscaledcov.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  write.table(combcov$allquants %>% dplyr::filter(gene == currgene) %>% 
                dplyr::select(transcript, method, TPM) %>%
                dplyr::mutate(TPM = round(TPM, 2)) %>% 
                tidyr::spread(method, TPM),
              file = paste0(outdir, "/tpm/", currgene, "_tpm.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  write.table(combcov$allquants %>% dplyr::filter(gene == currgene) %>% 
                dplyr::select(transcript, method, count) %>%
                dplyr::mutate(count = round(count, 2)) %>% 
                tidyr::spread(method, count),
              file = paste0(outdir, "/count/", currgene, "_count.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  saveRDS(jl, paste0(checkdir, "/", currgene, ".rds"))
}, mc.preschedule = FALSE, mc.cores = ncores)

saveRDS(NULL, paste0(checkdir, "/", basename(gene), ".rds"))
sessionInfo()
date()

