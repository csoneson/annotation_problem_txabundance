## ==================================================================================== ##
##              compare observed and predicted junction coverages                       ##
## ==================================================================================== ##
define obsvspredrule
figures/observed_vs_predicted_junction_coverage/observed_vs_predicted_junction_coverage_$(1)$(2).rds: output/$(1)$(2)_combined_coverages_with_scores.rds \
Rscripts/plot_observed_vs_predicted_junction_coverage.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$<' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,Salmon0.11,StringTie' uniqjuncreadsthreshold=$(uniqjuncreadsthreshold) fracuniqjuncreadsthreshold=0.75 outrds='$$@'" Rscripts/plot_observed_vs_predicted_junction_coverage.R Rout/plot_observed_vs_predicted_junction_coverage_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call obsvspredrule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call obsvspredrule,$(notdir $(F)),_stringtie_tx)))
$(foreach F,$(fastqfilesreal),$(eval $(call obsvspredrule,$(notdir $(F)),_chess)))

## ==================================================================================== ##
##           compare predicted coverage patterns between samples                        ##
## ==================================================================================== ##
figures/predicted_coverage_pattern_comparison/predicted_coverage_pattern_comparison_20151016.A-Cortex_RNA_20170918.A-WT_4.rds: \
alpine/20151016.A-Cortex_RNA/alpine_predicted_coverage.rds \
alpine/20170918.A-WT_4/alpine_predicted_coverage.rds Rscripts/plot_compare_predicted_coverage_patterns.R
	mkdir -p $(@D)
	$(R) "--args predcov1='$(word 1,$^)' predcov2='$(word 2,$^)' samplename1='Cortex' samplename2='HAP1' outrds='$@'" Rscripts/plot_compare_predicted_coverage_patterns.R Rout/plot_compare_predicted_coverage_patterns_20151016.A-Cortex_RNA_20170918.A-WT_4.Rout

## ==================================================================================== ##
##                            plot gene scores                                          ##
## ==================================================================================== ##
define plotscorerule
figures/gene_scores/gene_scores_$(1)$(2).rds: output/$(1)$(2)_combined_coverages_with_scores.rds \
Rscripts/plot_score_distribution.R Rscripts/define_plot_colors.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$(word 1,$$^)' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,Salmon0.11,StringTie' uniqjuncreadsthresholds='0,$(uniqjuncreadsthreshold)' outrds='$$@'" Rscripts/plot_score_distribution.R Rout/plot_score_distribution_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call plotscorerule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call plotscorerule,$(notdir $(F)),_stringtie_tx)))
$(foreach F,$(fastqfilesreal),$(eval $(call plotscorerule,$(notdir $(F)),_chess)))

## ==================================================================================== ##
##                   compare scores from different annotations                          ##
## ==================================================================================== ##
define scorecompannotchessrule
figures/comparison_scores_chess_ensembl/comparison_scores_chess_ensembl_$(1).rds: \
output/$(1)_chess_combined_coverages_with_scores.rds \
output/$(1)_combined_coverages_with_scores.rds \
reference_chess/chess2.0_assembly_fixed_tx2gene_withsymbol.rds \
Rscripts/plot_ensembl_vs_chess_scores.R
	mkdir -p $$(@D)
	$(R) "--args scorerdsensembl='output/$(1)_combined_coverages_with_scores.rds' scorerdschess='output/$(1)_chess_combined_coverages_with_scores.rds' convtablechess='reference_chess/chess2.0_assembly_fixed_tx2gene_withsymbol.rds' uniqjuncreadsthreshold=25 uniqjuncfracthreshold=0.75 outrds='$$@'" Rscripts/plot_ensembl_vs_chess_scores.R Rout/plot_ensembl_vs_chess_scores_$(1).Rout
endef
$(foreach F,$(fastqfilesreal),$(eval $(call scorecompannotchessrule,$(notdir $(F)))))

define scorecompannotstringtierule
figures/comparison_scores_stringtie_tx_ensembl/comparison_scores_stringtie_tx_ensembl_$(1).rds: \
output/$(1)_stringtie_tx_combined_coverages_with_scores.rds \
output/$(1)_combined_coverages_with_scores.rds \
reference/$(1)_stringtie_tx_tx2gene_withsymbol.rds \
reference/$(1)_stringtie_tx_tx2symbol.rds \
Rscripts/plot_ensembl_vs_stringtie_scores.R Rscripts/define_plot_colors.R
	mkdir -p $$(@D)
	$(R) "--args scorerdsensembl='output/$(1)_combined_coverages_with_scores.rds' scorerdsstringtie='output/$(1)_stringtie_tx_combined_coverages_with_scores.rds' convtablestringtie='reference/$(1)_stringtie_tx_tx2gene_withsymbol.rds' convtablestringtietx='reference/$(1)_stringtie_tx_tx2symbol.rds' uniqjuncreadsthreshold=25 uniqjuncfracthreshold=0.75 outrds='$$@'" Rscripts/plot_ensembl_vs_stringtie_scores.R Rout/plot_ensembl_vs_stringtie_scores_$(1).Rout
endef
$(foreach F,$(fastqfilesreal),$(eval $(call scorecompannotstringtierule,$(notdir $(F)))))

## ==================================================================================== ##
##                            genewise results/plots                                    ##
## ==================================================================================== ##
#define geneplotrule
#output_genewise/$(1)$(2)/check/$(3).rds: STARbigwig/$(1)$(2)_Aligned.sortedByCoord.out.bw output/$(1)$(2)_combined_coverages_with_scores.rds \
#$(gvizgenemodels) $(5) Rscripts/plot_genewise_results.R
#	mkdir -p output_genewise/$(1)$(2)/plots
#	mkdir -p output_genewise/$(1)$(2)/jcov
#	mkdir -p output_genewise/$(1)$(2)/tpm
#	mkdir -p output_genewise/$(1)$(2)/count
#	mkdir -p output_genewise/$(1)$(2)/check
#	$(R) "--args gene='$(3)' bigwig='$$(word 1,$$^)' bigwignanopore='$(5)' genemodels='$(gvizgenemodels)' scorerds='$$(word 2,$$^)' ncores=$(4) outdir='output_genewise/$(1)$(2)' libid='$(1)_' checkdir='output_genewise/$(1)$(2)/check'" Rscripts/plot_genewise_results.R Rout/plot_genewise_results_$(1)$(2)_$(3).Rout
#endef
#$(foreach G,$(genes_to_plot),$(foreach F,20151016.A-Cortex_RNA,$(eval $(call geneplotrule,$(F),,$(G),1,))))
#$(foreach G,$(genes_to_plot),$(foreach F,20170918.A-WT_4,$(eval $(call geneplotrule,$(F),,$(G),1,minimap2genomebigwig/20171207_1645_p2557_4017_2_ALLREADS.pass_minimap2_genome_s.bw))))
#$(foreach G,$(genes_to_plot),$(foreach F,sim_misannotated_utr_1,$(eval $(call geneplotrule,$(F),,$(G),1,))))

define geneplotrule2
figures/genewise_summary$(2)/$(1)$(2)_$(3).png: STARbigwig$(2)/$(1)_Aligned.sortedByCoord.out.bw output/$(1)$(2)_combined_coverages_with_scores.rds \
$(4) Rscripts/plot_genewise_summary.R
	mkdir -p $$(@D)
	$(R) "--args usegene='$(3)' bigwig='$$(word 1,$$^)' genemodels='$(4)' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,Salmon0.11,StringTie' scorerds='$$(word 2,$$^)' outpng='$$@'" Rscripts/plot_genewise_summary.R Rout/plot_genewise_summary_$(1)$(2)_$(3).Rout
endef
$(foreach G,$(genes_to_plot_summary),$(foreach F,$(fastqfilesreal),$(eval $(call geneplotrule2,$(notdir $(F)),,$(G),$(gvizgenemodels)))))
$(foreach G,$(genes_to_plot_summary_chess),$(foreach F,$(fastqfilesreal),$(eval $(call geneplotrule2,$(notdir $(F)),_chess,$(G),$(gvizgenemodels_chess)))))
$(foreach F,$(fastqfilesreal),$(foreach G,$(genes_to_plot_summary_stringtie$(notdir $(F))),$(eval $(call geneplotrule2,$(notdir $(F)),_stringtie_tx,$(G),stringtie/$(notdir $(F))/$(notdir $(F))_gviz_genemodels.rds))))

## ==================================================================================== ##
##                    correlation with inferential variance                             ##
## ==================================================================================== ##
define infvarrule
figures/correlation_with_inferential_variance/correlation_with_inferential_variance_$(1)$(2).rds: \
output/$(1)$(2)_combined_coverages_with_scores.rds $(3)/$(1)/quant.sf $(4) Rscripts/plot_correlation_with_inferential_variance.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$(word 1,$$^)' salmondir='$(3)/$(1)' uniqjuncreadsthreshold=$(uniqjuncreadsthreshold) tx2gene='$(4)' outrds='$$@'" Rscripts/plot_correlation_with_inferential_variance.R Rout/plot_correlation_with_inferential_variance_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call infvarrule,$(notdir $(F)),,salmon/cDNAncRNA,$(tx2gene))))

## ==================================================================================== ##
##                              correlation plots                                       ##
## ==================================================================================== ##
## Correlation between scores from different methods
define corrmethodrule
figures/correlation_between_methods/correlation_between_methods_$(1)$(2).rds: output/$(1)$(2)_combined_coverages_with_scores.rds \
Rscripts/plot_correlation_between_methods.R Rscripts/helper_plot_functions.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$(word 1,$$^)' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,Salmon0.11,StringTie' uniqjuncreadsthreshold=$(uniqjuncreadsthreshold) outrds='$$@'" Rscripts/plot_correlation_between_methods.R Rout/plot_correlation_between_methods_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call corrmethodrule,$(notdir $(F)),)))

## Correlation between scores from nanopore and illumina
define corrnanomethodrule
figures/correlation_between_nanopore_and_illumina_scores/correlation_between_nanopore_and_illumina_scores_$(1)$(2).rds: output/$(1)$(2)_combined_coverages_with_scores.rds Rscripts/plot_nanopore_vs_illumina_scores.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$(word 1,$$^)' targetmethod='WubMinimap2Nanopore' uniqjuncreadsthreshold=$(uniqjuncreadsthreshold) outrds='$$@'" Rscripts/plot_nanopore_vs_illumina_scores.R Rout/plot_nanopore_vs_illumina_scores_$(1)$(2).Rout
endef
$(eval $(call corrnanomethodrule,20170918.A-WT_4,))

## Correlation between estimated and true transcript abundances for simulated data
define corrtruthrule
figures/correlation_with_true_abundances/correlation_with_true_abundances_$(1)$(2).rds: output/$(1)$(2)_combined_coverages_with_scores.rds \
$(3) $(4) Rscripts/plot_estimated_abundance_accuracy.R Rscripts/helper_plot_functions.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$(word 1,$$^)' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,Salmon0.11,StringTie' truthrda='$(3)' truthmodgenesrds='$(4)' outrds='$$@'" Rscripts/plot_estimated_abundance_accuracy.R Rout/plot_estimated_abundance_accuracy_$(1)$(2).Rout
endef
$(eval $(call corrtruthrule,sim_misannotated_utr_1,,simulation/misannotated_utr/sim_counts_matrix.rda,simulation/misannotated_utr/sim_misannotated_utr_1_modified_genes.rds))
$(foreach F,$(fastqfilesreal),$(eval $(call corrtruthrule,$(notdir $(F)),,,)))

## Correlation between the within-gene Salmon/SalmonCDS abundance correlation and the Salmon score
define corrsalmonrule
figures/association_exoncdscorrelation_score/association_exoncdscorrelation_score_$(1)$(2).rds: output/$(1)$(2)_combined_coverages_with_scores.rds \
Rscripts/plot_association_exoncdscorrelation_score.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$(word 1,$$^)' outrds='$$@'" Rscripts/plot_association_exoncdscorrelation_score.R Rout/plot_association_exoncdscorrelation_score_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call corrsalmonrule,$(notdir $(F)),)))

## ==================================================================================== ##
##                         performance on simulated data                                ##
## ==================================================================================== ##
define simplotrule
figures/performance_simulated_data/performance_simulated_data_$(1)$(2).rds: output/$(1)$(2)_combined_coverages_with_scores.rds \
$(3) $(4) Rscripts/plot_performance_simdata.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$(word 1,$$^)' truthrda='$(3)' truthmodgenesrds='$(4)' gtf='$(5)' uniqjuncreadsthreshold=$(uniqjuncreadsthreshold) nthreads=$(nthreads) outrds='$$@'" Rscripts/plot_performance_simdata.R Rout/plot_performance_simdata_$(1)$(2).Rout
endef
$(eval $(call simplotrule,sim_misannotated_utr_1,,simulation/misannotated_utr/sim_counts_matrix.rda,simulation/misannotated_utr/sim_misannotated_utr_1_modified_genes.rds,$(gtf)))

## ==================================================================================== ##
##                       plot different annotation systems                              ##
## ==================================================================================== ##
define compareannotrule
figures/comparison_annotation_catalogs/annotation_comparison_$(1)_$(2).png: 
	mkdir -p $$(@D)
	$(R) "--args usegene='$(2)' bigwig='$(3)' baseannot='$(4)' basename='$(5)' annot1='$(6)' name1='$(7)' annot2='$(8)' name2='$(9)' outpng='$$@'" Rscripts/compare_annotation_catalogs.R Rout/compare_annotation_catalogs_$(1)_$(2)_$(5)_$(7)_$(9).Rout
endef
$(foreach F,$(fastqfilesreal),$(foreach G,$(genes_to_plot_summary),$(eval $(call compareannotrule,$(notdir $(F)),$(G),STARbigwig/$(notdir $(F))_Aligned.sortedByCoord.out.bw,$(gvizgenemodels),EnsemblGRCh38.90,$(gvizgenemodels_chess),CHESS,stringtie/$(notdir $(F))/$(notdir $(F))_gviz_genemodels.rds,StringTie))))

figures/ensembl_vs_chess_annotation_characteristics/ensembl_vs_chess_annotation_characteristics.rds: \
output/gene_characteristics.rds output/gene_characteristics_chess.rds \
reference_chess/chess2.0_assembly_fixed_tx2gene_withsymbol.rds Rscripts/plot_ensembl_vs_chess_annotation_characteristics.R
	mkdir -p $(@D)
	$(R) "--args genecharensembl='output/gene_characteristics.rds' genecharchess='output/gene_characteristics_chess.rds' convtablechess='reference_chess/chess2.0_assembly_fixed_tx2gene_withsymbol.rds' outrds='$@'" Rscripts/plot_ensembl_vs_chess_annotation_characteristics.R Rout/plot_ensembl_vs_chess_annotation_characteristics.Rout









