################################################################################
##                                                                            ##
## Plot observed vs predicted coverage for each junction                      ##
##                                                                            ##
## Inputs:                                                                    ##
## * scorerds: list containing abundance estimates and characteristics for    ##
##             junctions, transcripts and genes, as well as gene scores.      ##
##             Generated by calculate_gene_scores.R                           ##
## * quantmethods: string containing the quantification methods to consider,  ##
##                 separated by commas (no spaces)                            ##
## * uniqjuncreadsthreshold: the total number of uniquely mapping junction    ##
##                           reads (in a gene), will be used to stratify the  ##
##                           junctions                                        ##
## * fracuniqjuncreadsthreshold: the fraction of reads across a junction that ##
##                               map uniquely, will be used to stratify the   ##
##                               junctions                                    ##
## * outrds: output rds file. The name will be used to determine the name of  ##
##           the output figures.                                              ##
##                                                                            ##
## Outputs:                                                                   ##
## * A png figure illustrating the accuracy in "detection" (presence/absence  ##
##   calls) of junctions                                                      ##
## * A png figure showing the correlation between predicted and observed      ##
##   junction coverages                                                       ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

quantmethods <- strsplit(quantmethods, ",")[[1]]

print(scorerds)
print(quantmethods)
print(uniqjuncreadsthreshold)
print(fracuniqjuncreadsthreshold)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

## Read combined coverage table
x <- readRDS(scorerds)
x <- x$junctions

## Extract required information
y <- x %>% dplyr::select(gene, junctionid, method, pred.cov, scaled.cov, 
                         uniqreads, mmreads, fracunique) %>%
  tidyr::gather(covtype, coverage, pred.cov, scaled.cov) %>%
  dplyr::mutate(covtype = replace(covtype, covtype == "pred.cov", "Predicted coverage"),
                covtype = replace(covtype, covtype == "scaled.cov", "Scaled predicted coverage")) %>%
  dplyr::filter(method %in% quantmethods)

## Remove a few extreme outliers
print(sum(y$coverage >= (2 * max(y$uniqreads))))
y <- y %>% dplyr::filter(coverage < (2 * max(uniqreads)))

## ========================================================================== ##
## Compare predicted coverage to observed coverage in terms of being =0 or >0
z <- y %>% dplyr::filter(covtype == "Predicted coverage") %>%
  dplyr::mutate(coverage = round(coverage)) %>%
  dplyr::mutate(highfracunique = fracunique > fracuniqjuncreadsthreshold) %>%
  dplyr::group_by(gene, method) %>% 
  dplyr::mutate(totaluniqreads = sum(uniqreads)) %>% dplyr::ungroup() %>%
  dplyr::mutate(hightotaluniqreads = (totaluniqreads > uniqjuncreadsthreshold))

## Overall
g1 <- z %>% dplyr::group_by(method) %>% 
  dplyr::summarize(obs0pred0 = mean(coverage == 0 & uniqreads == 0),
                   obspospred0 = mean(coverage == 0 & uniqreads > 0),
                   obs0predlow = mean(coverage > 0 & coverage <= 5 & uniqreads == 0),
                   obs0predhigh = mean(coverage > 5 & uniqreads == 0),
                   obspospredpos = mean(coverage > 0 & uniqreads > 0)) %>%
  tidyr::gather(classif, fraction, -method) %>%
  dplyr::mutate(classif = replace(classif, classif == "obs0pred0", "Obs = 0, Pred = 0"),
                classif = replace(classif, classif == "obspospred0", "Obs > 0, Pred = 0"),
                classif = replace(classif, classif == "obs0predlow", "Obs = 0, 0 < Pred <= 5"),
                classif = replace(classif, classif == "obs0predhigh", "Obs = 0, Pred > 5"),
                classif = replace(classif, classif == "obspospredpos", "Obs > 0, Pred > 0")) %>%
  dplyr::mutate(classif = factor(classif, levels = c("Obs > 0, Pred = 0", "Obs = 0, Pred > 5",
                                                     "Obs = 0, 0 < Pred <= 5",
                                                     "Obs = 0, Pred = 0", "Obs > 0, Pred > 0"))) %>%
  ggplot(aes(x = method, y = fraction, fill = classif)) + 
  geom_bar(position = "fill", stat = "identity") + theme_bw() + 
  scale_fill_manual(values = c(`Obs > 0, Pred > 0` = "#000099", 
                               `Obs = 0, Pred = 0` = "#8080ff",
                               `Obs = 0, 0 < Pred <= 5` = "#ffc6b3",
                               `Obs = 0, Pred > 5` = "#ff6633",
                               `Obs > 0, Pred = 0` = "#992600"), name = "") + 
  xlab("") + ylab("Fraction of junctions") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## Split by total number of unique reads
g2 <- z %>% dplyr::group_by(method, hightotaluniqreads) %>% 
  dplyr::summarize(obs0pred0 = mean(coverage == 0 & uniqreads == 0),
                   obspospred0 = mean(coverage == 0 & uniqreads > 0),
                   obs0predlow = mean(coverage > 0 & coverage <= 5 & uniqreads == 0),
                   obs0predhigh = mean(coverage > 5 & uniqreads == 0),
                   obspospredpos = mean(coverage > 0 & uniqreads > 0)) %>%
  tidyr::gather(classif, fraction, -method, -hightotaluniqreads) %>%
  dplyr::mutate(classif = replace(classif, classif == "obs0pred0", "Obs = 0, Pred = 0"),
                classif = replace(classif, classif == "obspospred0", "Obs > 0, Pred = 0"),
                classif = replace(classif, classif == "obs0predlow", "Obs = 0, 0 < Pred <= 5"),
                classif = replace(classif, classif == "obs0predhigh", "Obs = 0, Pred > 5"),
                classif = replace(classif, classif == "obspospredpos", "Obs > 0, Pred > 0")) %>%
  dplyr::mutate(classif = factor(classif, levels = c("Obs > 0, Pred = 0", "Obs = 0, Pred > 5",
                                                     "Obs = 0, 0 < Pred <= 5",
                                                     "Obs = 0, Pred = 0", "Obs > 0, Pred > 0"))) %>%
  dplyr::mutate(hightotaluniqreads = replace(hightotaluniqreads, hightotaluniqreads == "TRUE", 
                                             paste0("Total unique junction reads > ", 
                                                    uniqjuncreadsthreshold)),
                hightotaluniqreads = replace(hightotaluniqreads, hightotaluniqreads == "FALSE", 
                                             paste0("Total unique junction reads <= ", 
                                                    uniqjuncreadsthreshold))) %>%
  ggplot(aes(x = method, y = fraction, fill = classif)) + 
  geom_bar(position = "fill", stat = "identity") + theme_bw() + 
  scale_fill_manual(values = c(`Obs > 0, Pred > 0` = "#000099", 
                               `Obs = 0, Pred = 0` = "#8080ff",
                               `Obs = 0, 0 < Pred <= 5` = "#ffc6b3",
                               `Obs = 0, Pred > 5` = "#ff6633",
                               `Obs > 0, Pred = 0` = "#992600"), name = "") + 
  xlab("") + ylab("Fraction of junctions") + facet_wrap(~ hightotaluniqreads, nrow = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## Split by fraction unique
g3 <- z %>% dplyr::group_by(method, highfracunique) %>% 
  dplyr::summarize(obs0pred0 = mean(coverage == 0 & uniqreads == 0),
                   obspospred0 = mean(coverage == 0 & uniqreads > 0),
                   obs0predlow = mean(coverage > 0 & coverage <= 5 & uniqreads == 0),
                   obs0predhigh = mean(coverage > 5 & uniqreads == 0),
                   obspospredpos = mean(coverage > 0 & uniqreads > 0)) %>%
  tidyr::gather(classif, fraction, -method, -highfracunique) %>%
  dplyr::mutate(classif = replace(classif, classif == "obs0pred0", "Obs = 0, Pred = 0"),
                classif = replace(classif, classif == "obspospred0", "Obs > 0, Pred = 0"),
                classif = replace(classif, classif == "obs0predlow", "Obs = 0, 0 < Pred <= 5"),
                classif = replace(classif, classif == "obs0predhigh", "Obs = 0, Pred > 5"),
                classif = replace(classif, classif == "obspospredpos", "Obs > 0, Pred > 0")) %>%
  dplyr::mutate(classif = factor(classif, levels = c("Obs > 0, Pred = 0", "Obs = 0, Pred > 5",
                                                     "Obs = 0, 0 < Pred <= 5",
                                                     "Obs = 0, Pred = 0", "Obs > 0, Pred > 0"))) %>%
  dplyr::mutate(highfracunique = replace(highfracunique, highfracunique == "TRUE", 
                                             paste0("Fraction unique reads > ",
                                                    fracuniqjuncreadsthreshold)),
                highfracunique = replace(highfracunique, highfracunique == "FALSE", 
                                             paste0("Fraction unique reads <= ",
                                                    fracuniqjuncreadsthreshold))) %>%
  ggplot(aes(x = method, y = fraction, fill = classif)) + 
  geom_bar(position = "fill", stat = "identity") + theme_bw() + 
  scale_fill_manual(values = c(`Obs > 0, Pred > 0` = "#000099", 
                               `Obs = 0, Pred = 0` = "#8080ff",
                               `Obs = 0, 0 < Pred <= 5` = "#ffc6b3",
                               `Obs = 0, Pred > 5` = "#ff6633",
                               `Obs > 0, Pred = 0` = "#992600"), name = "") + 
  xlab("") + ylab("Fraction of junctions") + facet_wrap(~ highfracunique, nrow = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

png(gsub("\\.rds$", "_detection.png", outrds), width = 7, height = 9, unit = "in", res = 300)
cowplot::plot_grid(g1, 
                   g2 + guides(fill = FALSE),
                   g3 + guides(fill = FALSE),
                   ncol = 1, labels = c("A", "B", "C"))
dev.off()

## ========================================================================== ##
## Scatter plots
## ========================================================================== ##
## Calculate correlations
corrs <- y %>% dplyr::group_by(method, covtype) %>% 
  dplyr::summarize(pearson = signif(cor(uniqreads, coverage, method = "pearson", 
                                        use = "pairwise.complete.obs"), 3),
                   spearman = signif(cor(uniqreads, coverage, method = "spearman", 
                                         use = "pairwise.complete.obs"), 3))

gg <- ggplot(y, aes(x = uniqreads, y = coverage)) + 
  geom_abline(intercept = 0, slope = 1, color = "black") + 
  geom_point(alpha = 0.3, size = 0.3, aes(color = (fracunique < fracuniqjuncreadsthreshold))) + 
  geom_label(data = corrs, x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1, 
             aes(label = paste0("Pearson: ", pearson, "\nSpearman: ", spearman))) + 
  facet_grid(method ~ covtype) + 
  xlab("Number of uniquely mapped reads spanning junction") + 
  ylab("Predicted number of reads spanning junction") + 
  scale_color_manual(name = paste0("Fraction\nuniquely mapping\nreads < ", 
                                   fracuniqjuncreadsthreshold), 
                     values = c(`TRUE` = "red", `FALSE` = "black")) + 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) + 
  theme_bw()

## Original scale
png(gsub("rds$", "png", outrds), width = 7, height = 15, unit = "in", res = 300)
print(gg)
dev.off()

## Square-root scale
png(gsub("\\.rds$", "_sqrt.png", outrds), width = 7, height = 15, unit = "in", res = 300)
print(gg + scale_x_sqrt() + scale_y_sqrt())
dev.off()

## MA-style
png(gsub("\\.rds$", "_sqrt_ma.png", outrds), width = 7, height = 15, unit = "in", res = 300)
print(ggplot(y, aes(x = (sqrt(uniqreads) + (sqrt(coverage)))/2, 
                    y = sqrt(coverage) - sqrt(uniqreads))) + 
        geom_abline(intercept = 0, slope = 0, color = "black") + 
        geom_point(alpha = 0.3, size = 0.3, 
                   aes(color = (fracunique < fracuniqjuncreadsthreshold))) + 
        facet_grid(method ~ covtype) + 
        xlab("(sqrt(Observed)+sqrt(Predicted))/2") + 
        ylab("sqrt(Predicted)-sqrt(Observed)") + 
        scale_color_manual(name = paste0("Fraction\nuniquely mapping\nreads < ", 
                                         fracuniqjuncreadsthreshold), 
                           values = c(`TRUE` = "red", `FALSE` = "black")) + 
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) + 
        theme_bw())
dev.off()

## hexbin
png(gsub("\\.rds$", "_sqrt_hex.png", outrds), width = 7, height = 15, unit = "in", res = 300)
print(ggplot(y, aes(x = uniqreads, y = coverage)) + 
        geom_abline(intercept = 0, slope = 1, color = "black") + 
        geom_hex(bins = 75) + 
        scale_fill_gradient(name = "", low = "blue", high = "red", trans = "log") + 
        geom_label(data = corrs, x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1, 
                   aes(label = paste0("Pearson: ", pearson, "\nSpearman: ", spearman))) + 
        facet_grid(method ~ covtype) + 
        xlab("Number of uniquely mapped reads spanning junction") + 
        ylab("Predicted number of reads spanning junction") + 
        theme_bw() + scale_x_sqrt() + scale_y_sqrt() + 
        theme(legend.position = "none"))
dev.off()

saveRDS(list(junctionscatter = gg), outrds)

date()
sessionInfo()
