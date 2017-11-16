args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gene)  ## gene of interest, or file listing collection of genes (one per row)
print(bigwig)  ## bigwig file for visualization
print(genemodels)  ## gene models etc (output of alpine_prepare_for_comparison.R)
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

## Determine which gene(s) to investigate
if (file.exists(gene)) {
  genes <- unlist(read.delim(gene, as.is = TRUE, header = FALSE))
} else {
  genes <- gene
}

replna <- function(x) {
  x[is.na(x)] <- 0
  x
}

## Investigate each gene
mclapply(genes, function(currgene) {
  jl <- genemodels$jcovscaled %>% dplyr::filter(gene == currgene) %>%
    dplyr::mutate(pred.cov = replna(pred.cov))
  
  jl0 <- jl %>% dplyr::select(seqnames, start, end, width, strand) %>%
    dplyr::distinct() %>% dplyr::arrange(start) %>% 
    dplyr::mutate(junctionid = paste0("J", seq_len(length(start))))
    
  jl <- jl %>% dplyr::left_join(genemodels$jcov, by = c("seqnames", "start", "end")) %>%
    dplyr::mutate(uniqreads = replace(uniqreads, is.na(uniqreads), 0),
                  mmreads = replace(mmreads, is.na(mmreads), 0)) %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(scaledcoverage = pred.cov/sum(pred.cov, na.rm = TRUE) * 
                    sum(uniqreads, na.rm = TRUE)) %>%
    dplyr::mutate(score = round(sum(abs(uniqreads - scaledcoverage), na.rm = TRUE)/
                                  sum(uniqreads, na.rm = TRUE), 2)) %>% 
    dplyr::mutate(methodscore = paste0(method, " (", score, ")")) %>%
    dplyr::left_join(jl0) %>% dplyr::ungroup() %>% 
    dplyr::select(junctionid, everything())

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
                tidyr::spread(methodscore, scaledcoverage) %%
                dplyr::arrange(start),
              file = paste0(outdir, "/jcov/", currgene, "_jcov.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  write.table(genemodels$allquants %>% dplyr::filter(gene == currgene) %>% 
                dplyr::select(transcript, method, TPM) %>%
                dplyr::mutate(TPM = round(TPM, 2)) %>% 
                tidyr::spread(method, TPM),
              file = paste0(outdir, "/tpm/", currgene, "_tpm.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  write.table(genemodels$allquants %>% dplyr::filter(gene == currgene) %>% 
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

