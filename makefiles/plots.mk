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

# ## ==================================================================================== ##
# ##                            select genes and plot                                     ##
# ## ==================================================================================== ##
# ## Generate text file with genes to investigate further
# define genelistrule
# gene_selection/$(1)/genes_to_run.txt: alpine/$(1)/alpine_predicted_coverage.rds $(tx2geneext) Rscripts/list_genes_to_run.R
# 	$(R) "--args inrds='$$<' tx2gene='$$(word 2,$$^)' outtxt='$$@'" Rscripts/list_genes_to_run.R Rout/list_genes_to_run_$(1).Rout
# endef
# $(foreach F,$(fastqfiles),$(eval $(call genelistrule,$(notdir $(F)))))

# ## TODO: FIX FOR STRINGTIE_TX
# ## Predict coverage and compare to observed junction coverage
# ## "gene" can be either a gene ID or a text file with a list of genes to investigate
# define alpinepredrule
# alpine_check/$(1)/$(notdir $(2)).rds: $(gvizgenemodels) \
# STARbigwig/$(1)_Aligned.sortedByCoord.out.bw alpine/$(1)/alpine_combined_coverages.rds \
# $(2) Rscripts/alpine_compare_coverage.R Rscripts/plot_tracks.R
# 	mkdir -p $$(@D)
# 	mkdir -p alpine_out/$(1)/plots
# 	mkdir -p alpine_out/$(1)/jcov
# 	mkdir -p alpine_out/$(1)/tpm
# 	mkdir -p alpine_out/$(1)/count	
# 	$(R) "--args gene='$(2)' bigwig='$$(word 2,$$^)' ncores=$(3) genemodels='$$(word 1,$$^)' combcovrds='$$(word 3,$$^)' outdir='alpine_out/$(1)' libid='$(1)_' checkdir='$$(@D)'" Rscripts/alpine_compare_coverage.R Rout/alpine_compare_coverage_$(1)_$(notdir $(2)).Rout
# endef
# $(foreach F,$(fastqfiles),$(eval $(call alpinepredrule,$(notdir $(F)),gene_selection/$(notdir $(F))/genes_to_run.txt,25)))
# $(foreach F,$(fastqfiles),$(eval $(call alpinepredrule,$(notdir $(F)),gene_selection/$(notdir $(F))/subset_genes_to_run.txt,25)))

# ## ==================================================================================== ##
# ##                            general summary plots                                     ##
# ## ==================================================================================== ##

# ## Deviation between predicted and observed coverages across data sets
# figures/deviation_predicted_observed_coverage_across_datasets.rds: alpine/20151016.A-Cortex_RNA/alpine_combined_coverages.rds \
# alpine/20170918.A-WT_4/alpine_combined_coverages.rds Rscripts/deviation_predicted_observed_coverage_across_datasets.R
# 	mkdir -p $(@D)
# 	$(R) "--args combcovrds1='$(word 1,$^)' combcovrds2='$(word 2,$^)' outrds='$@'" Rscripts/deviation_predicted_observed_coverage_across_datasets.R Rout/deviation_predicted_observed_coverage_across_datasets.Rout

# ## Correlation between scores from different methods
# define corrmethodrule
# figures/score_correlation_between_methods_$(1).rds: alpine/$(1)/alpine_combined_coverages.rds Rscripts/score_correlation_between_methods.R
# 	mkdir -p $$(@D)
# 	$(R) "--args combcovrds='$$(word 1,$$^)' outrds='$$@'" Rscripts/score_correlation_between_methods.R Rout/score_correlation_between_methods_$(1).Rout
# endef
# $(foreach F,$(fastqfiles),$(eval $(call corrmethodrule,$(notdir $(F)))))
