## ==================================================================================== ##
##              compare observed and predicted junction coverages                       ##
## ==================================================================================== ##
define obsvspredrule
figures/observed_vs_predicted_junction_coverage/observed_vs_predicted_junction_coverage_$(1)$(2).rds: output/$(1)$(2)_combined_coverages_with_scores.rds \
Rscripts/plot_observed_vs_predicted_junction_coverage.R
	mkdir -p $$(@D)
	$(R) "--args combcovrds='$$<' quantmethods='hera,kallisto,RSEM,Salmon,SalmonBWA,SalmonCDS,StringTie' outrds='$$@'" Rscripts/plot_observed_vs_predicted_junction_coverage.R Rout/plot_observed_vs_predicted_junction_coverage_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call obsvspredrule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call obsvspredrule,$(notdir $(F)),_stringtie_tx)))

## ==================================================================================== ##
##           compare predicted coverage patterns between samples                        ##
## ==================================================================================== ##
figures/predicted_coverage_pattern_comparison/predicted_coverage_pattern_comparison_20151016.A-Cortex_RNA_20170918.A-WT_4.rds: \
alpine/20151016.A-Cortex_RNA/alpine_predicted_coverage.rds \
alpine/20170918.A-WT_4/alpine_predicted_coverage.rds Rscripts/plot_compare_predicted_coverage_patterns.R
	mkdir -p $(@D)
	$(R) "--args predcov1='$(word 1,$^)' predcov2='$(word 2,$^)' outrds='$@'" Rscripts/plot_compare_predicted_coverage_patterns.R Rout/plot_compare_predicted_coverage_patterns_20151016.A-Cortex_RNA_20170918.A-WT_4.Rout

## ==================================================================================== ##
##                            plot gene scores                                          ##
## ==================================================================================== ##
define plotscorerule
figures/gene_scores/gene_scores_$(1)$(2).rds: output/$(1)$(2)_combined_coverages_with_scores.rds \
Rscripts/plot_score_distribution.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$(word 1,$$^)' quantmethods='hera,kallisto,RSEM,Salmon,SalmonBWA,SalmonCDS,StringTie' outrds='$$@'" Rscripts/plot_score_distribution.R Rout/plot_score_distribution_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call plotscorerule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call plotscorerule,$(notdir $(F)),_stringtie_tx)))

## ==================================================================================== ##
##                            genewise results/plots                                    ##
## ==================================================================================== ##
define geneplotrule
output_genewise/$(1)$(2)/check/$(3).rds: STARbigwig/$(1)$(2)_Aligned.sortedByCoord.out.bw output/$(1)$(2)_combined_coverages_with_scores.rds \
$(gvizgenemodels) Rscripts/plot_genewise_results.R
	mkdir -p output_genewise/$(1)$(2)/plots
	mkdir -p output_genewise/$(1)$(2)/jcov
	mkdir -p output_genewise/$(1)$(2)/tpm
	mkdir -p output_genewise/$(1)$(2)/count
	mkdir -p output_genewise/$(1)$(2)/check
	$(R) "--args gene='$(3)' bigwig='$$(word 1,$$^)' genemodels='$(gvizgenemodels)' combcovrds='$$(word 2,$$^)' ncores=$(4) outdir='output_genewise/$(1)$(2)' libid='$(1)_' checkdir='output_genewise/$(1)$(2)/check'" Rscripts/plot_genewise_results.R Rout/plot_genewise_results_$(1)$(2)_$(3).Rout
endef
$(foreach G,$(genes_to_plot),$(foreach F,$(fastqfiles),$(eval $(call geneplotrule,$(notdir $(F)),,$(G),1))))

## ==================================================================================== ##
##                    correlation with inferential variance                             ##
## ==================================================================================== ##
define infvarrule
figures/correlation_with_inferential_variance/correlation_with_inferential_variance_$(1)$(2).rds: \
output/$(1)$(2)_combined_coverages_with_scores.rds $(3)/$(1)/quant.sf $(4) Rscripts/plot_correlation_with_inferential_variance.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$(word 1,$$^)' salmondir='$(3)/$(1)' tx2gene='$(4)' outrds='$$@'" Rscripts/plot_correlation_with_inferential_variance.R Rout/plot_correlation_with_inferential_variance_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call infvarrule,$(notdir $(F)),,salmon/cDNAncRNA,$(tx2gene))))

## ==================================================================================== ##
##                            general summary plots                                     ##
## ==================================================================================== ##
## Correlation between scores from different methods
define corrmethodrule
figures/correlation_between_methods/correlation_between_methods_$(1)$(2).rds: output/$(1)$(2)_combined_coverages_with_scores.rds Rscripts/plot_correlation_between_methods.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$(word 1,$$^)' quantmethods='hera,kallisto,RSEM,Salmon,SalmonBWA,SalmonCDS,StringTie' outrds='$$@'" Rscripts/plot_correlation_between_methods.R Rout/plot_correlation_between_methods_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call corrmethodrule,$(notdir $(F)),)))
