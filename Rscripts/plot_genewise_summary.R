################################################################################
##                                                                            ##
## Generate reduced summary plots for one gene                                ##
##                                                                            ##
## Inputs:                                                                    ##
## * gene: gene of interest, or file listing genes of interest (one per row)  ##
## * bigwig: bigwig file of alignments for visualization                      ##
## * genemodels: gene models for Gviz                                         ##
## * quantmethods: string containing the quantification methods to consider,  ##
##                 separated by commas (no spaces)                            ##
## * scorerds: object with junction coverages, transcript abundances and gene ##
##             scores                                                         ##
## * outfile: output pdf file                                                 ##
##                                                                            ##
## Outputs:                                                                   ##
## * A png figure with coverage patterns, gene models and junction scores for ##
##   each gene of interest, as well as TPM estimates                          ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

quantmethods <- strsplit(quantmethods, ",")[[1]]

print(usegene)
print(bigwig)
print(genemodels) 
print(quantmethods)
print(scorerds)
print(outpng)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(rtracklayer)
  library(Gviz)
  library(scatterpie)
})

options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)

## Read gene models for Gviz plot (pregenerated from gtf to save time) and
## quantifications
genemodels <- readRDS(genemodels)$genemodels_exon
scores <- readRDS(scorerds)

jl <- scores$junctions %>% dplyr::filter(gene == usegene & 
                                           method %in% quantmethods)

bwfiles <- structure(bigwig, names = "")
bwcond <- structure("g1", names = "")

## Gene model track
gm <- subset(genemodels, tolower(gene) == tolower(usegene) | tolower(gene_name) == tolower(usegene))
gm <- subset(gm, gene == gene[1])  ## Select only one gene if there are many with the same name
id <- unique(gm$gene_name)
idshow <- paste0(id, " (", unique(gm$gene), ")")
show_chr <- unique(seqnames(gm))[1]
gm <- subset(gm, seqnames == show_chr)
min_coord <- min(start(gm)) - 0.2*(max(end(gm)) - min(start(gm)))
max_coord <- max(end(gm)) + 0.05*(max(end(gm)) - min(start(gm)))
gm$transcript <- factor(gm$transcript, levels = unique(gm$transcript))

## Other features in the considered region
gmo <- genemodels[overlapsAny(genemodels,
                              GRanges(seqnames = show_chr,
                                      ranges = IRanges(start = min_coord,
                                                       end = max_coord),
                                      strand = "*"))]
gmo <- gmo[!(gmo %in% gm)]
gmo <- reduce(gmo)

## Define colors
muted <- c("#DC050C","#E8601C","#7BAFDE","#1965B0","#B17BA6",
           "#882E72","#F1932D","#F6C141","#F7EE55","#4EB265",
           "#90C987","#CAEDAB","#777777")
txs <- levels(gm$transcript)
ncols <- nlevels(gm$transcript)
cols <- colorRampPalette(muted)(ncols)
names(cols) <- txs

grtr <- GeneRegionTrack(gm, showId = TRUE, col = NULL, fill = cols[gm$transcript],
                        name = "", col.title = "black", 
                        background.title = "transparent")
grtr2 <- GeneRegionTrack(gmo, showId = TRUE, col = "black", fill = "white",
                         name = "", col.title = "black", showId = FALSE,
                         background.title = "transparent")

gtr <- GenomeAxisTrack()

tracks <- c(gtr, grtr, grtr2)

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
tracks <- c(multiTracks_rnaseq, tracks)

rn <- round(1e6*runif(1))
png(paste0("gviz", rn, ".png"), width = 10.5, height = 5.25, unit = "in", res = 400)
plotTracks(tracks, chromosome = show_chr, 
           from = min_coord, to = max_coord, main = idshow, 
           min.width = 0, min.distance = 0, collapse = FALSE)
dev.off()

tpms <- ggplot(scores$transcripts %>% dplyr::filter(gene == usegene & method %in% quantmethods) %>%
                 dplyr::mutate(transcript = factor(transcript, levels = levels(gm$transcript))),
               aes(x = method, y = TPM, fill = transcript)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = cols, name = "") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
        legend.text = element_text(size = 7))

for (tt in txs) {
  jl[[tt]] <- grepl(tt, jl$transcript)
}
jl[, txs] <- sweep(jl[, txs], 1, rowSums(jl[, txs]), "/")

jcov <- ggplot() + geom_abline(intercept = 0, slope = 1) + 
  geom_scatterpie(aes(x = scaled.cov, y = uniqreads, r = max(scaled.cov)/13), 
                  cols = txs, data = jl, color = NA) + 
  facet_wrap(~ methodscore, nrow = 2) + 
  coord_equal(ratio = 1) + 
  expand_limits(x = range(c(jl$scaled.cov, jl$uniqreads)), y = range(c(jl$scaled.cov, jl$uniqreads))) + 
  scale_fill_manual(values = cols, name = "") + 
  xlab("Scaled predicted coverage") + 
  ylab("Number of uniquely mapped reads") + 
  theme_bw() + theme(strip.text = element_text(size = 7),
                     legend.text = element_text(size = 7))

# jcov <- ggplot(jl, aes(x = scaled.cov, y = uniqreads)) + 
#   geom_point(size = 3) + 
#   facet_wrap(~ methodscore, nrow = 2) + 
#   geom_abline(intercept = 0, slope = 1) + 
#   xlab("Scaled predicted coverage") + 
#   ylab("Number of uniquely mapped reads") + 
#   theme_bw() + theme(strip.text = element_text(size = 7))

png(outpng, width = 12, height = 10, unit = "in", res = 400)
print(cowplot::plot_grid(ggdraw() + draw_image(paste0("gviz", rn, ".png")),
                         plot_grid(tpms + theme(legend.position = "none"), 
                                   jcov + theme(legend.position = "none"), 
                                   nrow = 1, labels = c("B", "C"), rel_widths = c(0.9, 1)), 
                         get_legend(tpms + theme(legend.direction = "horizontal",
                                                 legend.justification = "center",
                                                 legend.box.just = "bottom") + 
                                      guides(fill = guide_legend(nrow = ifelse(length(txs) <= 8, 1, 2)))),
                         ncol = 1, rel_heights = c(1, 0.8, 0.1), 
                         labels = c("A", "", "")))
dev.off()

unlink(paste0("gviz", rn, ".png"))
date()
sessionInfo()