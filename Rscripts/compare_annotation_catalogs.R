################################################################################
##                                                                            ##
## Visualize multiple annotation catalogs, optionally together with a         ##
## coverag track                                                              ##
##                                                                            ##
## Inputs:                                                                    ##
## * usegene: gene of interest                                                ##
## * bigwig: bigwig file of alignments for visualization (can be '')          ##
## * baseannot: "base" gene models                                            ##
## * basename: track name for baseannot                                       ##
## * annot1,annot2: additional gene models                                    ##
## * name1,name2: track names for additional gene models                      ##
## * outpng: output png file                                                  ##
##                                                                            ##
## Outputs:                                                                   ##
## * A png figure with gene models                                            ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(usegene)
print(bigwig)
print(baseannot)
print(basename)
print(annot1)
print(name1)
print(annot2)
print(name2)
print(outpng)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(rtracklayer)
  library(Gviz)
})

options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)

## Read annotation files
baseannot <- readRDS(baseannot)$genemodels_exon
annot1 <- readRDS(annot1)$genemodels_exon
annot2 <- readRDS(annot2)$genemodels_exon

bwfiles <- structure(bigwig, names = "")
bwcond <- structure("g1", names = "")

gtr <- GenomeAxisTrack()

## Find region of interest (defined by gene in base annotation)
gm <- subset(baseannot, tolower(gene) == tolower(usegene) | tolower(gene_name) == tolower(usegene))
gm <- subset(gm, gene == gene[1])  ## Select only one gene if there are many with the same name
id <- unique(gm$gene_name)
idshow <- paste0(id, " (", unique(gm$gene), ")")
show_chr <- unique(seqnames(gm))[1]
gm <- subset(gm, seqnames == show_chr)
min_coord <- min(start(gm)) - 0.2*(max(end(gm)) - min(start(gm)))
max_coord <- max(end(gm)) + 0.05*(max(end(gm)) - min(start(gm)))

## Gene model tracks
grtrbase <- GeneRegionTrack(baseannot, showId = TRUE, col = NULL, fill = "#E8601C",
                            name = basename, col.title = "black", 
                            background.title = "white", min.height = 15)
grtr1 <- GeneRegionTrack(annot1, showId = TRUE, col = NULL, fill = "#7BAFDE",
                         name = name1, col.title = "black", 
                         background.title = "white", min.height = 15)
grtr2 <- GeneRegionTrack(annot2, showId = TRUE, col = NULL, fill = "#90C987",
                         name = name2, col.title = "black", 
                         background.title = "white", min.height = 15)

## Coverage track
multiTracks_rnaseq <- lapply(1:length(bwfiles), function(i) {
  assign(paste0("rnaseqtr", i), 
         DataTrack(range = bwfiles[i],
                   type = "histogram",
                   chromosome = unique(seqnames(gm)),
                   col.title = "black",
                   fill = "grey",
                   col = "grey",
                   col.histogram = "grey",
                   fill.histogram = "grey",
                   cex.title = 0,
         ))
})
tracks <- c(multiTracks_rnaseq, gtr, grtrbase, grtr1, grtr2)

png(outpng, width = 10, height = 9, unit = "in", res = 300)
plotTracks(tracks, chromosome = show_chr, 
           from = min_coord, to = max_coord, main = idshow, 
           min.width = 0, min.distance = 0, collapse = FALSE)
dev.off()

date()
sessionInfo()



