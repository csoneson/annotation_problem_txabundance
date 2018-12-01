################################################################################
##                                                                            ##
## Compare gene scores calculated with different functions (g)                ##
##                                                                            ##
## Inputs:                                                                    ##
## * combcovrds: object with junction coverage information for all methods    ##
##               (output from combine_scaled_coverages.R)                     ##
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

print(combcovrds)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(GGally)
  library(cowplot)
})

## Define weight functions
glist <- list(
  gthr0.75 = function(omega) {
    sapply(omega, function(o) {
      if (is.na(o) || o >= 0.75) 1
      else 0
    })
  },
  gsigmoid = function(omega) {
    sapply(omega, function(o) {
      if (is.na(o)) 1
      else 1/(1 + exp(-(25*(o - 0.7))))
    })
  },
  gpwlinear = function(omega) {
    sapply(omega, function(o) {
      if (is.na(o)) 1
      else if (o < 0.6) 0
      else 2.5*o - 1.5
    })
  },
  glinear = function(omega) {
    sapply(omega, function(o) {
      if (is.na(o)) 1
      else o
    })
  },
  gconstant = function(omega) {
    sapply(omega, function(o) {
      1
    })
  }
)

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

## Calculate score
junccov <- junccov %>%
  dplyr::group_by(gene, method) %>% 
  dplyr::mutate(score_gthr0.75_b1 = 
                  junction_score(uniqreads, mmreads, pred.cov, 
                                 g = glist$gthr0.75, beta = 1),
                score_gsigmoid_b1 = 
                  junction_score(uniqreads, mmreads, pred.cov,
                                 g = glist$gsigmoid, beta = 1),
                score_glinear_b1 = 
                  junction_score(uniqreads, mmreads, pred.cov,
                                 g = glist$glinear, beta = 1),
                score_gconstant_b1 = 
                  junction_score(uniqreads, mmreads, pred.cov,
                                 g = glist$gconstant, beta = 1),
                score_gpwlinear_b1 = 
                  junction_score(uniqreads, mmreads, pred.cov,
                                 g = glist$gpwlinear, beta = 1)) %>%
  dplyr::mutate(scaled.cov_gthr0.75_b1 = 
                  scaled_coverage(uniqreads, mmreads, pred.cov, 
                                  g = glist$gthr0.75, beta = 1),
                scaled.cov_gsigmoid_b1 = 
                  scaled_coverage(uniqreads, mmreads, pred.cov,
                                  g = glist$gsigmoid, beta = 1),
                scaled.cov_glinear_b1 = 
                  scaled_coverage(uniqreads, mmreads, pred.cov,
                                  g = glist$glinear, beta = 1),
                scaled.cov_gconstant_b1 = 
                  scaled_coverage(uniqreads, mmreads, pred.cov,
                                  g = glist$gconstant, beta = 1),
                scaled.cov_gpwlinear_b1 = 
                  scaled_coverage(uniqreads, mmreads, pred.cov,
                                  g = glist$gpwlinear, beta = 1)) %>%
  dplyr::ungroup()

## Add score to gene table
genecov <- combcov$gene
genecov <- dplyr::left_join(genecov, 
                            junccov %>% dplyr::select(gene, method, score_gthr0.75_b1,
                                                      score_gsigmoid_b1, score_glinear_b1,
                                                      score_gconstant_b1, score_gpwlinear_b1) %>%
                              dplyr::distinct(),
                            by = c("gene", "method"))

lowerfun <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.3) + 
    geom_point(alpha = 0.3, size = 0.5)
}

## Pairs plot of gene scores
gp <- ggpairs(genecov %>% dplyr::filter(method == "Salmon" & uniqjuncreads > 25) %>% 
                dplyr::select(score_gthr0.75_b1, score_gsigmoid_b1, score_gpwlinear_b1,
                              score_glinear_b1, score_gconstant_b1) %>%
                dplyr::rename(gthr0.75_b1 = score_gthr0.75_b1, 
                              gsigmoid_b1 = score_gsigmoid_b1,
                              gpwlinear_b1 = score_gpwlinear_b1,
                              glinear_b1 = score_glinear_b1,
                              gconstant_b1 = score_gconstant_b1),
          lower = list(continuous = lowerfun)) + 
  theme_bw() + xlab("JCC score") + ylab("JCC score")
gp

## Illustration of weight functions
x <- seq(0, 1, length.out = 1000)
df1 <- do.call(dplyr::bind_rows, list(
  data.frame(x = x, y = glist$gthr0.75(x), g = "gthr0.75_b1", stringsAsFactors = FALSE),
  data.frame(x = x, y = glist$gsigmoid(x), g = "gsigmoid_b1", stringsAsFactors = FALSE),
  data.frame(x = x, y = glist$gpwlinear(x), g = "gpwlinear_b1", stringsAsFactors = FALSE),
  data.frame(x = x, y = glist$glinear(x), g = "glinear_b1", stringsAsFactors = FALSE),
  data.frame(x = x, y = glist$gconstant(x), g = "gconstant_b1", stringsAsFactors = FALSE)
)) %>% dplyr::mutate(g = factor(g, levels = c("gthr0.75_b1", "gsigmoid_b1", "gpwlinear_b1",
                                              "glinear_b1", "gconstant_b1")))
gg <- ggplot(df1, aes(x = x, y = y)) + geom_path() + facet_wrap(~ g, ncol = 1) + 
  theme_bw() + xlab(expression(omega)) + ylab(expression(g(omega)))
gg

png(gsub("rds$", "png", outrds), width = 10, height = 7, unit = "in", res = 400)
cowplot::plot_grid(gg, ggmatrix_gtable(gp), nrow = 1, 
                   labels = c("A", "B"), rel_widths = c(0.3, 1))
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
