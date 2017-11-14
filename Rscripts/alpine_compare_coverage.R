args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gene)  ## gene of interest, or file listing collection of genes (one per row)
print(bam)  ## bam file
print(bigwig)  ## bigwig file for visualization
print(genemodels)  ## gene models etc (output of alpine_prepare_for_comparison.R)
print(biasmodels)  ## bias model object (output of alpine_fitbiasmodel.R)
print(ncores)  ## number of cores for parallel computations
print(outdir)  ## output directory
print(checkdir)  ## directory to write (empty) rds files (time stamps)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(alpine))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))

source("Rscripts/plot_tracks.R")

## Read gene models for Gviz plot (pregenerated from gtf to save time) and
## quantifications
genemodels <- readRDS(genemodels)

# ## Read bias model parameters and exon-by-transcript objects
biasmodels <- readRDS(biasmodels)
fitpar <- biasmodels$fitpar
ebt0 <- biasmodels$ebt0
txps <- biasmodels$txps

# ## Estimate average fragment length
avefraglength <- sum(fitpar$`1`$fraglen.density$x * fitpar$`1`$fraglen.density$y/
                       sum(fitpar$`1`$fraglen.density$y))

# ## Load bam file 
bam.files <- bam
names(bam.files) <- "1"

## Determine which gene(s) to investigate
if (file.exists(gene)) {
  genes <- unlist(read.delim(gene, as.is = TRUE, header = FALSE))
} else {
  genes <- gene
}

## Investigate each gene
mclapply(genes, function(currgene) {
  ## Get transcripts for gene of interest
  txlist <- names(subset(txps, gene_id == currgene))
  names(txlist) <- txlist
  
  if (length(txlist) > 0) {
    
    ## Predict coverage for each transcript
    pred.cov <- lapply(txlist, function(tx) {
      message(tx)
      ## Get transcript model
      txmod <- ebt0[[tx]]
      
      pc <- tryCatch({
        m <- predictCoverage(gene = txmod,
                             bam.files = bam.files,
                             fitpar = fitpar,
                             genome = Hsapiens,
                             model.names = "all")
        ## Scale predicted coverage to agree with estimated counts
        m$`1`$pred.cov$scaled <- list()
        for (nm in unique(genemodels$quants$method)) {
          m$`1`$pred.cov$scaled[[nm]] <- m$`1`$pred.cov$all/sum(m$`1`$pred.cov$all) * 
            genemodels$quants$count[genemodels$quants$transcript == tx & 
                                      genemodels$quants$method == nm] * avefraglength
        }
        m
      }, error = function(e) NULL)
      pc
    })
    
    junctionlist <- lapply(txlist, function(tx) {
      txmod <- sort(ebt0[[tx]])
      junctions <- GenomicRanges::setdiff(range(txmod), txmod)
      junctioncov <- list()
      if (all(strand(txmod) == "+")) {
        junctionpos <- cumsum(width(txmod))
        junctionpos <- junctionpos[-length(junctionpos)]
        for (nm in unique(genemodels$quants$method)) {
          junctioncov[[nm]] <- as.numeric(pred.cov[[tx]]$"1"$pred.cov$scaled[[nm]])[junctionpos]
        }
      } else if (all(strand(txmod) == "-")) {
        junctionpos <- cumsum(width(rev(txmod)))
        junctionpos <- junctionpos[-length(junctionpos)]
        for (nm in unique(genemodels$quants$method)) {
          junctioncov[[nm]] <- rev(as.numeric(pred.cov[[tx]]$"1"$pred.cov$scaled[[nm]])[junctionpos])
        }
      } else {
        stop("Unknown or mixed strand")
      }
      for (nm in unique(genemodels$quants$method)) {
        mcols(junctions)[, nm] <- junctioncov[[nm]]
      }
      junctions
    })
    
    replna <- function(x) {
      x[is.na(x)] <- 0
      x
    }
    
    jl <- do.call(rbind, lapply(junctionlist, function(w) {
      as.data.frame(w) %>%
        tidyr::gather(method, coverage, -seqnames, -start, -end, -width, -strand)
      })) %>% 
      dplyr::group_by(seqnames, start, end, width, strand, method) %>%
      dplyr::summarize(coverage = sum(coverage, na.rm = TRUE)) %>%
      dplyr::mutate(coverage = replna(coverage)) %>%
      dplyr::ungroup()

    jl0 <- jl %>% dplyr::select(seqnames, start, end, width, strand) %>%
      dplyr::distinct() %>% dplyr::arrange(start) %>% 
      dplyr::mutate(junctionid = paste0("J", seq_len(length(start))))
    
    jl <- dplyr::left_join(jl, genemodels$jcov) %>%
      dplyr::mutate(uniqreads = replace(uniqreads, is.na(uniqreads), 0),
                    mmreads = replace(mmreads, is.na(mmreads), 0)) %>%
      dplyr::group_by(method) %>%
      dplyr::mutate(scaledcoverage = coverage/sum(coverage, na.rm = TRUE) * 
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
    
    write.table(jl %>% dplyr::select(-score, -coverage, -method) %>%
                  dplyr::mutate(scaledcoverage = round(scaledcoverage, 2)) %>% 
                  tidyr::spread(methodscore, scaledcoverage),
                file = paste0(outdir, "/jcov/", currgene, "_jcov.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    write.table(genemodels$quants %>% dplyr::filter(transcript %in% txlist) %>% 
                  dplyr::select(transcript, method, TPM) %>%
                  dplyr::mutate(TPM = round(TPM, 2)) %>% 
                  tidyr::spread(method, TPM),
                file = paste0(outdir, "/tpm/", currgene, "_tpm.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    write.table(genemodels$quants %>% dplyr::filter(transcript %in% txlist) %>% 
                  dplyr::select(transcript, method, count) %>%
                  dplyr::mutate(count = round(count, 2)) %>% 
                  tidyr::spread(method, count),
                file = paste0(outdir, "/count/", currgene, "_count.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    saveRDS(jl, paste0(checkdir, "/", currgene, ".rds"))
  }
}, mc.preschedule = FALSE, mc.cores = ncores)

saveRDS(NULL, paste0(checkdir, "/", basename(gene), ".rds"))
sessionInfo()
date()

