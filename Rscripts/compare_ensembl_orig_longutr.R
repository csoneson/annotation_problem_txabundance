################################################################################
##                                                                            ##
## Compare gene scores with original Ensembl annotation and with annotation   ##
## extended with long 3'UTRs                                                  ##
##                                                                            ##
## Inputs:                                                                    ##
## * combcovrdsorig: object with junction coverage information for all        ##
##                   methods with original Ensembl annotation                 ##
## * combcovrdsorig: object with junction coverage information for all        ##
##                   methods with extended Ensembl annotation                 ##
## * uniqjuncreadsthr: threshold on the total number of uniquely mapping      ##
##                     junction reads                                         ##
## * uniqjuncfractionthr: threshold on the fraction of uniquely mapping       ## 
##                        junction reads                                      ##
## * outrds: output file                                                      ##
##                                                                            ##
## Outputs:                                                                   ##
## * Plots comparing the JCC scores                                           ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(combcovrdsorig)
print(combcovrdslongutr)
print(uniqjuncreadsthr)
print(uniqjuncfractionthr)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

combcovorig <- readRDS(combcovrdsorig)
combcovlongutr <- readRDS(combcovrdslongutr)

orig_genes <- combcovorig$genes
longutr_genes <- combcovlongutr$genes

df_gene <- dplyr::full_join(
  orig_genes %>% dplyr::filter(uniqjuncreads > uniqjuncreadsthr & 
                                 uniqjuncfraction > uniqjuncfractionthr) %>% 
    dplyr::select(gene, method, score) %>%
    dplyr::rename(orig = score),
  longutr_genes %>% dplyr::filter(uniqjuncreads > uniqjuncreadsthr & 
                                    uniqjuncfraction > uniqjuncfractionthr) %>% 
    dplyr::select(gene, method, score) %>%
    dplyr::rename(longutr = score)
) %>% dplyr::filter(!is.na(orig) & !is.na(longutr)) %>%
  dplyr::mutate(longutr_minus_orig = longutr - orig) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(max_longutr_minus_orig = max(longutr_minus_orig)) %>%
  dplyr::ungroup()

png(gsub("rds", "png", outrds), width = 11, height = 4, unit = "in", res = 400)
ggplot(df_gene, aes(x = orig, y = longutr)) +
  geom_abline(slope = 1, intercept = 0) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~ method, nrow = 1) + theme_bw() + 
  xlab("Original annotation catalog") + 
  ylab("Annotation catalog extended with long 3'UTRs")
dev.off()

lapply(split(df_gene, f = df_gene$method), function(x) {
  x %>% dplyr::arrange(longutr - orig) %>%
    head(n = 5)
})

write.table(
  df_gene %>% dplyr::arrange(max_longutr_minus_orig),
  file = gsub("\\.rds", "_consistently_improved_in_longutr.txt", outrds),
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)

saveRDS(df_gene, file = outrds)
date()
sessionInfo()