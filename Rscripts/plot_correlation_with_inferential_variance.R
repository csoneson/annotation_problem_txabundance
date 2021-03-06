################################################################################
##                                                                            ##
## Plot scores against inferential variance, estimated via bootstrapping      ##
##                                                                            ##
## Inputs:                                                                    ##
## * scorerds: list containing abundance estimates and characteristics for    ##
##             junctions, transcripts and genes, as well as gene scores.      ##
##             Generated by calculate_gene_scores.R                           ##
## * salmondir: the path to a directory containing output from Salmon,        ##
##              including bootstrap estimates                                 ##
## * uniqjuncreadsthreshold: the total number of uniquely mapping junction    ##
##                           reads (in a gene), only genes with more than     ##
##                           this number will be used for the comparison      ##
## * tx2gene: a data frame mapping transcript IDs to gene IDs                 ##
## * outrds: output rds file. The name will be used to determine the name of  ##
##           the output figures.                                              ##
##                                                                            ##
## Outputs:                                                                   ##
## * A png figure comparing the score to the inferential coefficient of       ##
##   variation, for genes and transcripts.                                    ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(scorerds)
print(salmondir)
print(uniqjuncreadsthreshold)
print(tx2gene)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(tximport)
})

## Read scores
scores <- readRDS(scorerds)$genes

## Read bootstrap counts on transcript level
tx2gene <- readRDS(tx2gene)
salmonres <- tximport(paste0(salmondir, "/quant.sf"), type = "salmon", txOut = TRUE)
bootcounts <- salmonres$infReps[[1]]
rownames(bootcounts) <- rownames(salmonres$counts)

## Summarize counts on gene level
bootcounts_gene <- as.data.frame(bootcounts) %>% tibble::rownames_to_column("tx") %>%
  dplyr::left_join(tx2gene %>% dplyr::select(tx, gene), by = "tx") %>%
  dplyr::select(-tx) %>% dplyr::group_by(gene) %>% dplyr::summarize_all(funs(sum)) %>%
  as.data.frame() %>% tibble::column_to_rownames(var = "gene")
idx <- grep("^STRG\\.|^CHS\\.", rownames(bootcounts_gene), invert = TRUE)
rownames(bootcounts_gene)[idx] <- gsub("\\.[0-9]+$", "", rownames(bootcounts_gene)[idx])

## Calculate coefficient of variation for each transcript and gene
cv_tx <- apply(bootcounts, 1, sd)/apply(bootcounts, 1, mean)
cv_tx[is.na(cv_tx)] <- 0
cv_gene <- apply(bootcounts_gene, 1, sd)/apply(bootcounts_gene, 1, mean)
cv_gene[is.na(cv_gene)] <- 0

## Merge with score information. For transcripts, we use the score for the
## corresponding gene
df_tx <- data.frame(tx = names(cv_tx), CV = cv_tx, stringsAsFactors = FALSE) %>%
  dplyr::left_join(tx2gene %>% dplyr::select(tx, gene))

idx <- grep("^STRG\\.|^CHS\\.", df_tx$tx, invert = TRUE)
df_tx$tx[idx] <- gsub("\\.[0-9]+$", "", df_tx$tx[idx])
idx <- grep("^STRG\\.|^CHS\\.", df_tx$gene, invert = TRUE)
df_tx$gene[idx] <- gsub("\\.[0-9]+$", "", df_tx$gene[idx])

df_tx <- df_tx %>% 
  dplyr::inner_join(scores %>% dplyr::filter(method == "Salmon" & 
                                               uniqjuncreads >= uniqjuncreadsthreshold) %>%
                      dplyr::select(gene, score))

df_gene <- scores %>% dplyr::filter(method == "Salmon" & 
                                      uniqjuncreads >= uniqjuncreadsthreshold) %>% 
  dplyr::select(gene, score) %>%
  dplyr::inner_join(data.frame(gene = names(cv_gene), CV = cv_gene, stringsAsFactors = FALSE))

## Plot
png(gsub("rds$", "png", outrds), width = 6, height = 12, unit = "in", res = 300)
cowplot::plot_grid(
  ggplot(df_tx, aes(x = CV, y = score)) + 
    geom_point(alpha = 0.3, size = 1) + geom_smooth() + 
    theme_bw() + ggtitle("Transcript") + ylab("JCC score") + 
    theme(axis.title = element_text(size = 14)),
  
  ggplot(df_gene %>% dplyr::mutate(CV = replace(CV, CV > 1, 1)), aes(x = CV, y = score)) + 
    geom_point(alpha = 0.3, size = 1) + geom_smooth(data = df_gene %>% dplyr::filter(CV < 1)) + 
    theme_bw() + ggtitle("Gene") + ylab("JCC score") + 
    theme(axis.title = element_text(size = 14)),
  
  ncol = 1
)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
