################################################################################
##                                                                            ##
## Generate gene2tx file for RSEM                                             ##
##                                                                            ##
## Inputs:                                                                    ##
## * tx2gene: data frame with transcript-to-gene conversion information       ##
## * rsemgene2tx: output file name                                            ##
##                                                                            ##
## Outputs:                                                                   ##
## * File with gene-to-transcript conversion information for RSEM             ##
##                                                                            ##
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(tx2gene)
print(rsemgene2tx)

x <- readRDS(tx2gene)

write.table(x[, c("gene", "tx")], file = rsemgene2tx, 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

sessionInfo()
date()
