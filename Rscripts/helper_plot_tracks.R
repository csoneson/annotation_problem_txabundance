suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(GenomicRanges))

create_genemodels <- function(gtf_file, seltype = "exon") {
  genemodels <- import(gtf_file)
  if (all(c("transcript_id", "gene_id", "exon_id") %in% colnames(mcols(genemodels)))) {
    idx <- match(c("transcript_id", "gene_id", "exon_id"), colnames(mcols(genemodels)))
  } else {
    idx <- match(c("transcript_id", "gene_id", "exon_number"), colnames(mcols(genemodels)))
  }
  colnames(mcols(genemodels))[idx] <- c("transcript", "gene", "exon")
  mcols(genemodels)$symbol <- mcols(genemodels)$transcript
  genemodels <- subset(genemodels, type == seltype)
  
  genemodels
}

plot_tracks <- function(mygene = NULL, genemodels = NULL, genemodels2 = NULL, 
                        gtf_file = NULL, rnaseq_datafiles = NULL, 
                        rnaseq_condition = NULL, show_chr = NULL, 
                        min_coord = NULL, max_coord = NULL, 
                        pdf_filename = NULL, pdf_width = 7, pdf_height = 7, ...) {
  options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)

  if (is.null(genemodels) & is.null(gtf_file)) 
    stop("Either genemodels or gtf_file must be provided")
  if (is.null(genemodels))
    genemodels <- create_genemodels(gtf_file)
  if (!is.null(genemodels))
    stopifnot(all(c("transcript", "gene", "exon") %in% colnames(mcols(genemodels))))
  
  ## Create gene region tracks
  if (!is.null(mygene)) {
    ## Track 1
    gm <- subset(genemodels, tolower(gene) == tolower(mygene) | tolower(gene_name) == tolower(mygene))
    gm <- subset(gm, gene == gene[1])  ## Select only one gene if there are many with the same name
    id <- unique(gm$gene_name)
    idshow <- paste0(id, " (", unique(gm$gene), ")")
    show_chr <- unique(seqnames(gm))[1]
    gm <- subset(gm, seqnames == show_chr)
    min_coord <- min(start(gm)) - 0.15*(max(end(gm)) - min(start(gm)))
    max_coord <- max(end(gm)) + 0.05*(max(end(gm)) - min(start(gm)))
    
    ## Other features in the considered region
    gmo <- genemodels[overlapsAny(genemodels,
                                  GRanges(seqnames = show_chr,
                                          ranges = IRanges(start = min_coord,
                                                           end = max_coord),
                                          strand = "*"))]
    gmo <- gmo[!(gmo %in% gm)]
    
    ## Track 2
    if (!is.null(genemodels2)) {
      gm2 <- subset(genemodels2, tolower(gene) == tolower(mygene) | tolower(gene_name) == tolower(mygene))
      gm2 <- subset(gm2, gene == gm$gene[1])  ## Select only one gene if there are many with the same name
      gm2 <- subset(gm2, seqnames == show_chr)
      min_coord2 <- min(start(gm2)) - 0.15*(max(end(gm2)) - min(start(gm2)))
      max_coord2 <- max(end(gm2)) + 0.05*(max(end(gm2)) - min(start(gm2)))
      min_coord <- min(min_coord, min_coord2)
      max_coord <- max(max_coord, max_coord2)
      
      # gm2 <- genemodels2[overlapsAny(genemodels2, 
      #                                GRanges(seqnames = show_chr,
      #                                        ranges = IRanges(start = min_coord,
      #                                                         end = max_coord),
      #                                        strand = "*")), ]
      # ## Excluded features from other genes will be part of the "other genes" track
      # gm2other <- subset(gm2, gene != gm$gene[1])
      # mcols(gm2other) <- mcols(gm2other)[, match(colnames(mcols(gmo)), colnames(mcols(gm2other)))]
      # gmo <- c(gmo, gm2other)
      # ## Keep only excluded features from the current gene in this track
      # gm2 <- subset(gm2, gene == gm$gene[1])
    } else {
      gm2 <- GRanges()
    }
  } else {
    if (any(is.null(c(show_chr, min_coord, max_coord)))) {
      stop("Either mygene or genomic region must be provided")
    } else {
      id <- NULL
      idshow <- ""
      gm <- genemodels[overlapsAny(genemodels,
                                   GRanges(seqnames = show_chr, 
                                           ranges = IRanges(start = min_coord,
                                                            end = max_coord),
                                           strand = "*")), ]
      gm2 <- genemodels2[overlapsAny(genemodels2,
                                     GRanges(seqnames = show_chr, 
                                             ranges = IRanges(start = min_coord,
                                                              end = max_coord),
                                             strand = "*")), ]
      gmo <- NULL
    }
  }
  grtr <- GeneRegionTrack(gm, showId = TRUE, col = NULL, fill = "blue",
                          name = ifelse(!is.null(id), paste0(id, ", exons"), "Genes"), 
                          col.title = "black")
  grtr2 <- GeneRegionTrack(gmo, showId = TRUE, col = NULL, fill = "green",
                           name = "", 
                           col.title = "black")
  grtr3 <- GeneRegionTrack(gm2, showId = TRUE, col = NULL, fill = "orange",
                           name = ifelse(!is.null(id), paste0(id, ", CDSs"), "Genes"), 
                           col.title = "black")
  
  ## Create genome axis track
  gtr <- GenomeAxisTrack()
  
  if (length(gm) > 0) {
    tracks <- c(gtr, grtr, grtr3, grtr2)
  } else {
    tracks <- c(gtr)
  }
  
  threecols <- c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D", "black")
  # twocols <- c(rgb(11, 102, 254, maxColorValue = 255), 
  #              rgb(250, 0, 255, maxColorValue = 255))
  
  ## RNAseq tracks
  if (!is.null(rnaseq_datafiles)) {
    if (is.null(names(rnaseq_datafiles))) 
      stop("RNAseq input file list must be named")
    multiTracks_rnaseq <- lapply(1:length(rnaseq_datafiles), function(i) {
      assign(paste0("rnaseqtr", i), 
             DataTrack(range = rnaseq_datafiles[i],
                       type = "histogram",
                       name = names(rnaseq_datafiles)[i],
                       chromosome = unique(seqnames(gm)),
                       col.title = "black",
                       fill = threecols[(as.numeric(as.factor(rnaseq_condition)))[i]],
                       col = threecols[(as.numeric(as.factor(rnaseq_condition)))[i]],
                       col.histogram = threecols[(as.numeric(as.factor(rnaseq_condition)))[i]],
                       fill.histogram = threecols[(as.numeric(as.factor(rnaseq_condition)))[i]]
             ))
    })
    tracks <- c(multiTracks_rnaseq, tracks)
  }
  
  ## Plot tracks
  if (!is.null(pdf_filename))
    pdf(pdf_filename, width = pdf_width, height = pdf_height)
  
  plotTracks(tracks, chromosome = show_chr, 
             from = min_coord, to = max_coord, 
             main = ifelse(!is.null(id), idshow, ""), 
             min.width = 0, min.distance = 0, collapse = FALSE, ...)
  
  if (!is.null(pdf_filename))
    dev.off()
}