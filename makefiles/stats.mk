define covpredsumrule
stats/alpine_coverage_prediction_summary_$(1)$(2).txt: alpine/$(1)$(2)/alpine_predicted_coverage.rds \
Rscripts/summarize_predicted_coverage_alpine.R
	mkdir -p $$(@D)
	$(R) "--args covrds='$$<' outtxt='$$@'" Rscripts/summarize_predicted_coverage_alpine.R Rout/summarize_predicted_coverage_alpine_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call covpredsumrule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call covpredsumrule,$(notdir $(F)),_stringtie_tx)))
$(foreach F,$(fastqfiles),$(eval $(call covpredsumrule,$(notdir $(F)),_chess)))

## Extract genes with high scores
define highscorerule
stats/genes_with_high_score_$(1)$(2).txt: output/$(1)$(2)_combined_coverages_with_scores.rds \
Rscripts/get_highscore_genes.R
	mkdir -p $$(@D)
	$(R) "--args scorerds='$$<' quantmethods='hera,kallisto,RSEM,Salmon,SalmonSTAR,SalmonCDS,SalmonKeepDup,StringTie' uniqjuncreadsthreshold=$(uniqjuncreadsthreshold) outtxt='$$@'" Rscripts/get_highscore_genes.R Rout/get_highscore_genes_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call highscorerule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call highscorerule,$(notdir $(F)),_stringtie_tx)))
$(foreach F,$(fastqfiles),$(eval $(call highscorerule,$(notdir $(F)),_chess)))

## Characterize annotation catalogs (gtf files)
define annotcharrule
stats/annotation_characteristics$(1).txt: $(2) $(3) $(4) Rscripts/characterize_annotation.R
	mkdir -p $$(@D)
	$(R) "--args gtf='$(2)' txome='$(3)' idconvrds='$(4)' outtxt='$$@'" Rscripts/characterize_annotation.R Rout/characterize_annotation$(1).Rout
endef
$(eval $(call annotcharrule,_ensembl.38.90,$(gtf),$(txome),))
$(eval $(call annotcharrule,_chess2.0_assembly_fixed,$(gtf_chess),$(txome_chess),$(tx2gene_chess_withsymbol)))
$(foreach F,$(fastqfilesreal),$(eval $(call annotcharrule,_$(notdir $(F))_stringtie_tx,stringtie/$(notdir $(F))/$(notdir $(F))_filtered.gtf,stringtie/$(notdir $(F))/$(notdir $(F))_stringtie_tx.fa,reference/$(notdir $(F))_stringtie_tx_tx2gene_withsymbol.rds)))
