################################################################################
##                                                                            ##
## Count number of duplicated transcripts that are removed by Salmon          ##
##                                                                            ##
## Inputs:                                                                    ##
## * salmonindexdir: Salmon index directory                                   ##
## * tx2gene: A data frame mapping transcripts to genes                       ##
## * outtxt: Output text file                                                 ##
##                                                                            ##
## Outputs:                                                                   ##
## * A text file indicating the number of removed transcripts                 ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(salmonindexdir)
print(tx2gene)
print(outtxt)

x <- read.delim(paste0(salmonindexdir, "/duplicate_clusters.tsv"), 
                header = TRUE, as.is = TRUE)
tx2gene <- readRDS(tx2gene)
stopifnot(all(x$DuplicateTxp %in% tx2gene$tx))

sink(outtxt)
print(paste0("Number of removed transcripts: ", length(unique(x$DuplicateTxp))))
print(paste0("Number of genes corresponding to removed transcripts: ",
             length(unique(tx2gene$gene[tx2gene$tx %in% x$DuplicateTxp]))))

print(paste0("Total number of genes from the start: ", length(unique(tx2gene$gene))))
print(paste0("Total number of genes completely removed: ", length(unique(tx2gene$gene)) - length(unique(tx2gene$gene[!(tx2gene$tx %in% x$DuplicateTxp)]))))
sink()

date()
sessionInfo()

