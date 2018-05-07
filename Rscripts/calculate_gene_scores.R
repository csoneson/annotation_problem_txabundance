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

## Read combined coverage file
combcov <- readRDS(combcovrds)
combcov <- combcov$jcovscaled

## Calculate score. Consider only junctions with not too many multimapping reads
combcov <- combcov %>%
  dplyr::mutate(include.junction = mmreads/(uniqreads + mmreads) < mmfracthreshold) %>%
  dplyr::mutate(include.junction = replace(include.junction, is.na(include.junction), 1)) %>%
  dplyr::group_by(gene, method) %>% 
  dplyr::mutate(score = round(sum(abs(uniqreads - scaled.cov), na.rm = TRUE)/
                                sum(uniqreads, na.rm = TRUE), 2)) %>% 
  dplyr::mutate(methodscore = paste0(method, " (", score, ")")) %>%
  dplyr::mutate(scoreMM = round(sum((abs(uniqreads - scaled.cov.mm)) * include.junction, na.rm = TRUE)/
                                  sum(uniqreads * include.junction, na.rm = TRUE), 2)) %>% 
  dplyr::mutate(methodscoreMM = paste0(method, " (", scoreMM, ")")) %>%
  dplyr::ungroup()

saveRDS(combcov, file = outrds)

date()
sessionInfo()