## ==================================================================================== ##
##              compare observed and predicted junction coverages                       ##
## ==================================================================================== ##
define obsvspredrule
figures/observed_vs_predicted_junction_coverage/observed_vs_predicted_junction_coverage_$(1)$(2).rds: output/$(1)$(2)_combined_coverages_with_scores.rds \
Rscripts/plot_observed_vs_predicted_junction_coverage.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$<' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,StringTie' uniqjuncreadsthreshold=$(uniqjuncreadsthreshold) fracuniqjuncreadsthreshold=$(fracuniqjuncreadsthreshold) outrds='$$@'" Rscripts/plot_observed_vs_predicted_junction_coverage.R Rout/plot_observed_vs_predicted_junction_coverage_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call obsvspredrule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call obsvspredrule,$(notdir $(F)),_stringtie_tx)))
$(foreach F,$(fastqfilesreal),$(eval $(call obsvspredrule,$(notdir $(F)),_chess)))

define obsvspredpermrule
figures/observed_vs_predicted_junction_coverage/observed_vs_predicted_junction_coverage_$(1)$(2)_permuted.rds: output/$(1)$(2)_combined_coverages_permuted_with_scores.rds \
Rscripts/plot_observed_vs_predicted_junction_coverage.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$<' quantmethods='heraPermuted,kallistoPermuted,RSEMPermuted,SalmonPermuted,SalmonSTARPermuted,SalmonCDSPermuted,SalmonKeepDupPermuted,StringTiePermuted' uniqjuncreadsthreshold=$(uniqjuncreadsthreshold) fracuniqjuncreadsthreshold=$(fracuniqjuncreadsthreshold) outrds='$$@'" Rscripts/plot_observed_vs_predicted_junction_coverage.R Rout/plot_observed_vs_predicted_junction_coverage_$(1)$(2)_permuted.Rout
endef
$(eval $(call obsvspredpermrule,20170918.A-WT_4,))

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
	$(R) "--args scorerds='$$(word 1,$$^)' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,StringTie' uniqjuncreadsthresholds='0,$(uniqjuncreadsthreshold)' outrds='$$@'" Rscripts/plot_score_distribution.R Rout/plot_score_distribution_$(1)$(2).Rout
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
	$(R) "--args scorerdsensembl='output/$(1)_combined_coverages_with_scores.rds' scorerdschess='output/$(1)_chess_combined_coverages_with_scores.rds' convtablechess='reference_chess/chess2.0_assembly_fixed_tx2gene_withsymbol.rds' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,StringTie' uniqjuncreadsthreshold=$(uniqjuncreadsthreshold) uniqjuncfracthreshold=0.75 outrds='$$@'" Rscripts/plot_ensembl_vs_chess_scores.R Rout/plot_ensembl_vs_chess_scores_$(1).Rout
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
	$(R) "--args scorerdsensembl='output/$(1)_combined_coverages_with_scores.rds' scorerdsstringtie='output/$(1)_stringtie_tx_combined_coverages_with_scores.rds' convtablestringtie='reference/$(1)_stringtie_tx_tx2gene_withsymbol.rds' convtablestringtietx='reference/$(1)_stringtie_tx_tx2symbol.rds' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,StringTie' uniqjuncreadsthreshold=$(uniqjuncreadsthreshold) uniqjuncfracthreshold=0.75 outrds='$$@'" Rscripts/plot_ensembl_vs_stringtie_scores.R Rout/plot_ensembl_vs_stringtie_scores_$(1).Rout
endef
$(foreach F,$(fastqfilesreal),$(eval $(call scorecompannotstringtierule,$(notdir $(F)))))

## ==================================================================================== ##
##                            genewise results/plots                                    ##
## ==================================================================================== ##
define geneplotrule2
figures/genewise_summary$(2)/$(1)$(2)_$(3).png: STARbigwig$(5)/$(1)_Aligned.sortedByCoord.out.bw output/$(1)$(2)_combined_coverages_with_scores.rds \
$(4) Rscripts/plot_genewise_summary.R
	mkdir -p $$(@D)
	$(R) "--args usegene='$(3)' bigwig='$$(word 1,$$^)' genemodels='$(4)' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,StringTie' scorerds='$$(word 2,$$^)' outpng='$$@'" Rscripts/plot_genewise_summary.R Rout/plot_genewise_summary_$(1)$(2)_$(3).Rout
endef
$(foreach G,$(genes_to_plot_summary),$(foreach F,$(fastqfilesreal),$(eval $(call geneplotrule2,$(notdir $(F)),,$(G),$(gvizgenemodels),))))
$(foreach G,$(genes_to_plot_summary_chess),$(foreach F,$(fastqfilesreal),$(eval $(call geneplotrule2,$(notdir $(F)),_chess,$(G),$(gvizgenemodels_chess),_chess))))
$(foreach F,$(fastqfilesreal),$(foreach G,$(genes_to_plot_summary_stringtie$(notdir $(F))),$(eval $(call geneplotrule2,$(notdir $(F)),_stringtie_tx,$(G),stringtie/$(notdir $(F))/$(notdir $(F))_gviz_genemodels.rds,_stringtie_tx))))
$(foreach G,$(genes_to_plot_summary_longutr_added),$(foreach F,$(fastqfilesreal),$(eval $(call geneplotrule2,$(notdir $(F)),_longUTR_added,$(G),$(gvizgenemodels_longutr_added),))))

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
	$(R) "--args scorerds='$$(word 1,$$^)' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,StringTie' uniqjuncreadsthreshold=$(uniqjuncreadsthreshold) outrds='$$@'" Rscripts/plot_correlation_between_methods.R Rout/plot_correlation_between_methods_$(1)$(2).Rout
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
	$(R) "--args scorerds='$$(word 1,$$^)' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,StringTie' truthrda='$(3)' truthmodgenesrds='$(4)' outrds='$$@'" Rscripts/plot_estimated_abundance_accuracy.R Rout/plot_estimated_abundance_accuracy_$(1)$(2).Rout
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

## ==================================================================================== ##
##                compare scores with different weight functions                        ##
## ==================================================================================== ##
define compareweightrule
figures/comparison_scores_diff_weights/comparison_score_diff_weights_$(1)$(2).rds: output/$(1)$(2)_combined_coverages.rds \
Rscripts/compare_score_definitions.R
	mkdir -p $$(@D)
	$(R) "--args combcovrds='output/$(1)$(2)_combined_coverages.rds' outrds='$$@'" Rscripts/compare_score_definitions.R Rout/compare_score_definitions_$(1)$(2).Rout
endef
$(foreach F,$(fastqfilesreal),$(eval $(call compareweightrule,$(notdir $(F)),)))

## ==================================================================================== ##
##           plot derived junction coverages from STAR and HISAT2 alignments            ##
## ==================================================================================== ##
define starhisatjuncrule
figures/comparison_junccov_star_hisat2/comparison_junccov_star_hisat2_$(1).rds: \
STAR/$(1)/$(1)_Aligned.sortedByCoord.out.bam HISAT2/$(1)/$(1).bam Rscripts/plot_junc_cov_star_hisat2.R
	mkdir -p $$(@D)
	$(R) "--args junctioncovSTAR='STAR/$(1)/$(1)_SJ.out.tab' bamFileSTAR='STAR/$(1)/$(1)_Aligned.sortedByCoord.out.bam' bamFileHISAT2='HISAT2/$(1)/$(1).bam' outrds='$$@'" Rscripts/plot_junc_cov_star_hisat2.R Rout/plot_junc_cov_star_hisat2_$(1).Rout
endef
$(foreach F,$(fastqfilesreal),$(eval $(call starhisatjuncrule,$(notdir $(F)))))

## ==================================================================================== ##
##         compare scores between Ensembl annotation and Ensembl + long 3'UTRs          ##
## ==================================================================================== ##
define ensembllongutrscorerule
figures/comparison_orig_longutr_scores/comparison_orig_longutr_scores_$(1).rds: \
output/$(1)_combined_coverages_with_scores.rds output/$(1)_longUTR_added_combined_coverages_with_scores.rds \
Rscripts/compare_ensembl_orig_longutr.R
	mkdir -p $$(@D)
	$(R) "--args combcovrdsorig='$$(word 1,$$^)' combcovrdslongutr='$$(word 2,$$^)' uniqjuncreadsthr=$(uniqjuncreadsthreshold) uniqjuncfractionthr=$(fracuniqjuncreadsthreshold) outrds='$$@'" Rscripts/compare_ensembl_orig_longutr.R Rout/compare_ensembl_orig_longutr_$(1).Rout
endef
$(foreach F,$(fastqfilesreal), $(eval $(call ensembllongutrscorerule,$(notdir $(F)))))



