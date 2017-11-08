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

calc_prop_p <- function(coverage, uniqreads, tot_coverage, tot_reads) {
  sapply(seq_along(coverage), function(i) {
    prop.test(x = c(coverage[i], uniqreads[i]), n = c(tot_coverage[i], tot_reads[i]),
              alternative = "two.sided")$p.value
  })
}

## Read gene models for Gviz plot (pregenerated from gtf to save time) and
## quantifications
genemodels <- readRDS(genemodels)

## Read bias model parameters and exon-by-transcript objects
biasmodels <- readRDS(biasmodels)
fitpar <- biasmodels$fitpar
ebt0 <- biasmodels$ebt0
txps <- biasmodels$txps

## Estimate average fragment length
avefraglength <- sum(fitpar$`1`$fraglen.density$x * fitpar$`1`$fraglen.density$y/
                       sum(fitpar$`1`$fraglen.density$y))

## Load bam file 
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
        ## Scale predicted coverage to agree with Salmon's estimated count
        m$`1`$pred.cov$scaled <- list()
        for (nm in names(genemodels$quants)) {
          m$`1`$pred.cov$scaled[[nm]] <- m$`1`$pred.cov$all/sum(m$`1`$pred.cov$all) * 
            genemodels$quants[[nm]]$count[genemodels$quants[[nm]]$transcript == tx] * avefraglength
        }
        # m$`1`$pred.cov$all <- m$`1`$pred.cov$all/sum(m$`1`$pred.cov$all) * 
        #   genemodels$quantsf$NumReads[genemodels$quantsf$Name == tx] * avefraglength
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
        for (nm in names(genemodels$quants)) {
          junctioncov[[nm]] <- as.numeric(pred.cov[[tx]]$"1"$pred.cov$scaled[[nm]])[junctionpos]
        }
        # junctioncov <- as.numeric(pred.cov[[tx]]$"1"$pred.cov$all)[junctionpos]
      } else if (all(strand(txmod) == "-")) {
        junctionpos <- cumsum(width(rev(txmod)))
        junctionpos <- junctionpos[-length(junctionpos)]
        for (nm in names(genemodels$quants)) {
          junctioncov[[nm]] <- rev(as.numeric(pred.cov[[tx]]$"1"$pred.cov$scaled[[nm]])[junctionpos])
        }
        # junctioncov <- as.numeric(pred.cov[[tx]]$"1"$pred.cov$all)[junctionpos]
        # junctioncov <- rev(junctioncov)
      } else {
        stop("Unknown or mixed strand")
      }
      for (nm in names(genemodels$quants)) {
        mcols(junctions)[, paste0(nm, "_coverage")] <- junctioncov[[nm]]
      }
      # mcols(junctions)$coverage <- junctioncov
      junctions
    })
    
    replna <- function(x) {
      x[is.na(x)] <- 0
      x
    }
    
    covcols <- grep("_coverage", colnames(mcols(junctionlist[[1]])), value = TRUE)
    jl <- do.call(rbind, lapply(junctionlist, as.data.frame)) %>% 
      dplyr::group_by(seqnames, start, end, width, strand) %>%
      dplyr::summarize_at(.vars = covcols, .funs = sum, na.rm = TRUE) %>%
      # dplyr::summarize(coverage = sum(coverage, na.rm = TRUE)) %>% 
      ungroup() %>%
      dplyr::mutate_at(.vars = covcols, .funs = funs(replna))
      # dplyr::mutate(coverage = replace(coverage, is.na(coverage), 0))
    
    jl <- dplyr::left_join(jl, genemodels$jcov) %>%
      dplyr::mutate(uniqreads = replace(uniqreads, is.na(uniqreads), 0),
                    mmreads = replace(mmreads, is.na(mmreads), 0)) %>%
      dplyr::mutate_at(.vars = covcols, .funs = funs(./sum(., na.rm = TRUE) * sum(uniqreads, na.rm = TRUE))) %>%
      # dplyr::mutate(scaledcoverage = coverage/sum(coverage, na.rm = TRUE) * sum(uniqreads, na.rm = TRUE)) %>%
      # dplyr::mutate(tot_coverage = sum(coverage, na.rm = TRUE),
      #               tot_reads = sum(uniqreads, na.rm = TRUE)) %>%
      # dplyr::mutate(prop_pval = calc_prop_p(coverage, uniqreads, tot_coverage, tot_reads)) %>%
      dplyr::mutate(junctionid = paste0("J", seq_len(length(uniqreads)))) %>% 
      # dplyr::mutate(difference = uniqreads - scaledcoverage) %>%
      # dplyr::mutate(ranking = order(order(difference))) %>%
      dplyr::select(junctionid, everything())
    
    pdf(paste0(outdir, "/", currgene, ".pdf"), width = 12, height = 10)
    tryCatch({
      plot_tracks(mygene = currgene, genemodels = genemodels$genemodels_exon, 
                  genemodels2 = genemodels$genemodels_cds, 
                  gtf_file = NULL, rnaseq_datafiles = structure(bigwig, names = "s1"), 
                  rnaseq_condition = structure("g1", names = "s1"), show_chr = NULL, 
                  min_coord = NULL, max_coord = NULL, 
                  pdf_filename = NULL, pdf_width = 7, pdf_height = 7)
    }, error = function(e) message(e))
    
    grid.newpage()
    quants0 <- genemodels$quants
    for (nm in names(quants0)) {
      idx <- which(colnames(quants0[[nm]]) != "transcript")
      colnames(quants0[[nm]])[idx] <- paste0(nm, "_", colnames(quants0[[nm]])[idx])
    }
    quants0 <- lapply(quants0, function(x) x %>% dplyr::filter(transcript %in% txlist))
    quants0 <- Reduce(function(...) dplyr::full_join(..., by = "transcript"), quants0)
    grid.table(quants0 %>% dplyr::select(-StringTie_EffectiveLength))
    # grid.table(genemodels$quantsf %>% dplyr::filter(Name %in% txlist))
    
    grid.newpage()
    grid.table(jl %>% dplyr::select(junctionid, seqnames, start, end, width, strand, 
                                    uniqreads, mmreads, Salmon_coverage, kallisto_coverage,
                                    RSEM_coverage, StringTie_coverage))

    methods <- names(genemodels$quants)
    names(methods) <- methods
    genescore <- lapply(methods, function(nm) {
      round(sum(abs(jl$uniqreads - jl[, paste0(nm, "_coverage")]), na.rm = TRUE)/
              sum(jl$uniqreads, na.rm = TRUE), 2)
    })    
    # genescore <- round(sum(abs(jl$uniqreads - jl$scaledcoverage), 
    #                        na.rm = TRUE)/sum(jl$uniqreads, na.rm = TRUE), 2)
    
    plots <- lapply(methods, function(nm) {
      ggplot(jl, aes_string(x = paste0(nm, "_coverage"), y = "uniqreads", 
                            label = "junctionid")) + 
        geom_point(size = 4) + geom_label_repel() + 
        geom_abline(intercept = 0, slope = 1) + 
        ggtitle(paste0(nm, ", score = ", genescore[[nm]])) + 
        xlab("Scaled predicted coverage") + ylab("Number of uniquely mapped reads")
    })
    print(plot_grid(plotlist = plots))

    # print(ggplot(jl, aes(x = scaledcoverage, y = uniqreads, label = junctionid)) + 
    #         geom_point(size = 4) + geom_label_repel() + 
    #         geom_abline(intercept = 0, slope = 1) + 
    #         ggtitle(paste0("score = ", genescore)) + 
    #         xlab("Scaled predicted coverage") + ylab("Number of uniquely mapped reads"))
    dev.off()
    
    write.table(jl %>%
                  dplyr::mutate_at(.vars = covcols, .funs = funs(round(., 2))),
                file = paste0(outdir, "/", currgene, ".txt"),
                quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    # print(as.data.frame(jl))
    
    saveRDS(list(score = genescore), paste0(checkdir, "/", currgene, ".rds"))
  }
}, mc.preschedule = FALSE, mc.cores = ncores)

saveRDS(NULL, paste0(checkdir, "/", basename(gene), ".rds"))
sessionInfo()
date()

