## ==================================================================================== ##
##                            characterize genes                                        ##
## ==================================================================================== ##
## From annotation
output/characterize_genes.rds: $(gtf) $(txome) Rscripts/characterize_genes.R
	$(R) "--args gtf='$(gtf)' txome='$(txome)' outrds='$@'" Rscripts/characterize_genes.R Rout/characterize_genes.Rout

## Summarize gene expression from all methods
define combgexrule
alpine/$(1)$(2)/alpine_gene_expression.rds: alpine/$(1)$(2)/alpine_combined_coverages.rds \
Rscripts/combine_gene_expression_estimates.R
	$(R) "--args combcovrds='$$(word 1,$$^)' outrds='$$@'" Rscripts/combine_gene_expression_estimates.R Rout/combine_gene_expression_estimates_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call combgexrule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call combgexrule,$(notdir $(F)),_stringtie_tx)))

## ==================================================================================== ##
##                            plot gene scores                                          ##
## ==================================================================================== ##
define plotscorerule
alpine/$(1)$(2)/$(1)$(2)_gene_scores.rds: alpine/$(1)$(2)/alpine_gene_expression.rds \
output/characterize_genes.rds alpine/$(1)$(2)/alpine_combined_coverages.rds \
featureCounts/$(1)/$(1)_STAR_exons.txt \
featureCounts/$(1)/$(1)_STAR_introns.txt \
Rscripts/plot_score_distribution.R
	$(R) "--args covrds='$$(word 3,$$^)' gexrds='$$(word 1,$$^)' geneinfords='$$(word 2,$$^)' exoncountstxt='$$(word 4,$$^)' introncountstxt='$$(word 5,$$^)' outrds='$$@'" Rscripts/plot_score_distribution.R Rout/plot_score_distribution_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call plotscorerule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call plotscorerule,$(notdir $(F)),_stringtie_tx)))

## ==================================================================================== ##
##                            select genes and plot                                     ##
## ==================================================================================== ##
## Generate text file with genes to investigate further
define genelistrule
gene_selection/$(1)/genes_to_run.txt: alpine/$(1)/alpine_predicted_coverage.rds $(tx2geneext) Rscripts/list_genes_to_run.R
	$(R) "--args inrds='$$<' tx2gene='$$(word 2,$$^)' outtxt='$$@'" Rscripts/list_genes_to_run.R Rout/list_genes_to_run_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call genelistrule,$(notdir $(F)))))

## TODO: FIX FOR STRINGTIE_TX
## Predict coverage and compare to observed junction coverage
## "gene" can be either a gene ID or a text file with a list of genes to investigate
define alpinepredrule
alpine_check/$(1)/$(notdir $(2)).rds: $(gvizgenemodels) \
STARbigwig/$(1)_Aligned.sortedByCoord.out.bw alpine/$(1)/alpine_combined_coverages.rds \
$(2) Rscripts/alpine_compare_coverage.R Rscripts/plot_tracks.R
	mkdir -p $$(@D)
	mkdir -p alpine_out/$(1)/plots
	mkdir -p alpine_out/$(1)/jcov
	mkdir -p alpine_out/$(1)/tpm
	mkdir -p alpine_out/$(1)/count	
	$(R) "--args gene='$(2)' bigwig='$$(word 2,$$^)' ncores=$(3) genemodels='$$(word 1,$$^)' combcovrds='$$(word 3,$$^)' outdir='alpine_out/$(1)' libid='$(1)_' checkdir='$$(@D)'" Rscripts/alpine_compare_coverage.R Rout/alpine_compare_coverage_$(1)_$(notdir $(2)).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call alpinepredrule,$(notdir $(F)),gene_selection/$(notdir $(F))/genes_to_run.txt,25)))
$(foreach F,$(fastqfiles),$(eval $(call alpinepredrule,$(notdir $(F)),gene_selection/$(notdir $(F))/subset_genes_to_run.txt,25)))

## ==================================================================================== ##
##                            general summary plots                                     ##
## ==================================================================================== ##
## Compare predicted coverage patterns between samples
figures/compare_predicted_coverage_patterns_20151016.A-Cortex_RNA_20170918.A-WT_4.rds: alpine/20151016.A-Cortex_RNA/alpine_predicted_coverage.rds \
alpine/20170918.A-WT_4/alpine_predicted_coverage.rds Rscripts/compare_predicted_coverage_patterns.R
	mkdir -p $(@D)
	$(R) "--args predcov1='$(word 1,$^)' predcov2='$(word 2,$^)' outrds='$@'" Rscripts/compare_predicted_coverage_patterns.R Rout/compare_predicted_coverage_patterns_20151016.A-Cortex_RNA_20170918.A-WT_4.Rout

## Overall correlation between predicted and observed coverages
define corrobspredrule
figures/correlation_predicted_observed_coverage_$(1).rds: alpine/$(1)/alpine_combined_coverages.rds Rscripts/correlation_predicted_observed_coverage.R
	mkdir -p $$(@D)
	$(R) "--args combcovrds='$$(word 1,$$^)' outrds='$$@'" Rscripts/correlation_predicted_observed_coverage.R Rout/correlation_predicted_observed_coverage_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call corrobspredrule,$(notdir $(F)))))

## Deviation between predicted and observed coverages across data sets
figures/deviation_predicted_observed_coverage_across_datasets.rds: alpine/20151016.A-Cortex_RNA/alpine_combined_coverages.rds \
alpine/20170918.A-WT_4/alpine_combined_coverages.rds Rscripts/deviation_predicted_observed_coverage_across_datasets.R
	mkdir -p $(@D)
	$(R) "--args combcovrds1='$(word 1,$^)' combcovrds2='$(word 2,$^)' outrds='$@'" Rscripts/deviation_predicted_observed_coverage_across_datasets.R Rout/deviation_predicted_observed_coverage_across_datasets.Rout

