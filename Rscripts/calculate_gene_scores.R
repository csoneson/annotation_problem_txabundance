args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(combcovrds)
print(mmfracthreshold)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
})

## Define weight functions
gthr <- function(omega, thr = mmfracthreshold) {
  sapply(omega, function(o) {
    if (is.na(o) || o >= (1 - thr)) 1
    else 0
  })
}

## Define help function for calculating score
junction_score <- function(uniqreads, mmreads, predcovs, g, beta = 1) {
  omega <- uniqreads/(uniqreads + mmreads)
  omega[mmreads == 0] <- 1  ## if there are no multi-mapping reads, all reads are considered to be unique
  w1 <- (sum(g(omega) * uniqreads)/sum(g(omega) * predcovs)) ^ beta
  ## w1 can be non-numeric if all g(omega)=0 (not enough uniquely mapping reads
  ## for any junction) or if g(omega)*pred.cov=0 for all junctions, even if
  ## g(omega)!=0 for some of them (if the predicted coverage is 0 for a junction
  ## that has non-zero uniquely mapping reads). In both these cases, we don't
  ## scale the predicted coverage (i.e., we set w1=1).
  w1[is.na(w1)] <- 1
  w1[!is.finite(w1)] <- 1
  signif(sum(abs(w1 * g(omega) * predcovs - g(omega) * uniqreads))/sum(g(omega) * uniqreads), 2)
}

## Define corresponding help function for calculating scaled coverages
scaled_coverage <- function(uniqreads, mmreads, predcovs, g, beta = 1) {
  omega <- uniqreads/(uniqreads + mmreads)
  omega[mmreads == 0] <- 1  ## if there are no multi-mapping reads, all reads are considered to be unique
  w1 <- (sum(g(omega) * uniqreads)/sum(g(omega) * predcovs)) ^ beta
  ## w1 can be non-numeric if all g(omega)=0 (not enough uniquely mapping reads
  ## for any junction) or if g(omega)*pred.cov=0 for all junctions, even if
  ## g(omega)!=0 for some of them (if the predicted coverage is 0 for a junction
  ## that has non-zero uniquely mapping reads). In both these cases, we don't
  ## scale the predicted coverage (i.e., we set w1=1).
  w1[is.na(w1)] <- 1
  w1[!is.finite(w1)] <- 1
  w1 * predcovs
}

## Read combined coverage file
combcov <- readRDS(combcovrds)
junccov <- combcov$junctions

## Calculate score. Consider only junctions with not too many multimapping reads
junccov <- junccov %>%
  dplyr::mutate(fracunique = uniqreads/(uniqreads + mmreads)) %>% 
  dplyr::mutate(fracunique = replace(fracunique, is.na(fracunique), 1)) %>%
  dplyr::group_by(gene, method) %>% 
  dplyr::mutate(score = junction_score(uniqreads, mmreads, pred.cov, g = gthr, beta = 1)) %>%
  dplyr::mutate(scaled.cov = scaled_coverage(uniqreads, mmreads, pred.cov, g = gthr, beta = 1)) %>%
  dplyr::mutate(methodscore = paste0(method, " (", score, ")")) %>%
  dplyr::ungroup()

## Add score to gene table
genecov <- combcov$gene
genecov <- dplyr::left_join(genecov, 
                            junccov %>% dplyr::select(gene, method, score) %>%
                              dplyr::distinct(),
                            by = c("gene", "method"))

saveRDS(list(junctions = junccov, transcripts = combcov$transcripts,
             genes = genecov), file = outrds)

date()
sessionInfo()