## Fit bias model
define fitbiasrule
alpine/$(1)$(2)/alpine_fitbiasmodel.rds: $(3) STAR$(2)/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
Rscripts/alpine_fitbiasmodel.R
	mkdir -p $$(@D)
	$(R) "--args gtf='$(3)' bam='$$(word 2,$$^)' readlength=$(4) minsize=$(5) maxsize=$(6) organism='$(7)' genomeVersion='$(8)' version=$(9) outdir='$$(@D)'" Rscripts/alpine_fitbiasmodel.R Rout/alpine_fitbiasmodel_$(1)$(2).Rout
endef
$(eval $(call fitbiasrule,20151016.A-Cortex_RNA,,$(gtf),126,100,300,Homo_sapiens,GRCh38,90))
$(eval $(call fitbiasrule,20170918.A-WT_4,,$(gtf),151,140,450,Homo_sapiens,GRCh38,90))
$(eval $(call fitbiasrule,20151016.A-Cortex_RNA,_stringtie_tx,stringtie/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA_filtered.gtf,126,100,300,Homo_sapiens,GRCh38,90))
$(eval $(call fitbiasrule,20170918.A-WT_4,_stringtie_tx,stringtie/20170918.A-WT_4/20170918.A-WT_4_filtered.gtf,151,140,450,Homo_sapiens,GRCh38,90))

## Predict transcript and junction coverage profiles for all transcripts that have at least one 
## junction and are longer than the fragment length
define predcovrule
alpine/$(1)$(2)/alpine_predicted_coverage.rds: alpine/$(1)$(2)/alpine_fitbiasmodel.rds \
STAR$(2)/$(1)/$(1)_Aligned.sortedByCoord.out.bam Rscripts/alpine_get_predicted_coverage.R
	mkdir -p $$(@D)
	$(R) "--args bam='$$(word 2,$$^)' biasmodels='$$(word 1,$$^)' ncores=$(3) outrds='$$@'" Rscripts/alpine_get_predicted_coverage.R Rout/alpine_get_predicted_coverage_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call predcovrule,$(notdir $(F)),,25)))
$(foreach F,$(fastqfiles),$(eval $(call predcovrule,$(notdir $(F)),_stringtie_tx,25)))

## Scale junction coverage by transcript abundance estimates for each method
define juncscalerule
alpine/$(1)$(2)/scaled_junction_coverage_$(3).rds: alpine/$(1)$(2)/alpine_predicted_coverage.rds \
$(4) $(5) $(6) Rscripts/alpine_scale_junction_coverage.R
	mkdir -p $$(@D)
	$(R) "--args predcovrds='$$(word 1,$$^)' txquants='$(4)' quantreadscript='$(5)' tx2gene='$(6)' strandspec='$(7)' method='$(3)' outrds='$$@'" Rscripts/alpine_scale_junction_coverage.R Rout/alpine_scale_junction_coverage_$(1)$(2)_$(3).Rout
endef
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,Salmon,salmon/cDNAncRNA/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,SalmonBWA,salmonbwa/cDNAncRNA/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,SalmonCDS,salmon/cds/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,kallisto,kallisto/cDNAncRNA/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_kallisto.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,RSEM,RSEM/cDNAncRNA/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA.isoforms.results,Rscripts/read_quant_rsem.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,hera,hera/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_hera.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,StringTie,stringtie_onlyref/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA.gtf,Rscripts/read_quant_stringtie.R,$(tx2geneext),yes))

$(eval $(call juncscalerule,20170918.A-WT_4,,Salmon,salmon/cDNAncRNA/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20170918.A-WT_4,,SalmonBWA,salmonbwa/cDNAncRNA/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20170918.A-WT_4,,SalmonCDS,salmon/cds/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20170918.A-WT_4,,kallisto,kallisto/cDNAncRNA/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_kallisto.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20170918.A-WT_4,,RSEM,RSEM/cDNAncRNA/20170918.A-WT_4/20170918.A-WT_4.isoforms.results,Rscripts/read_quant_rsem.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20170918.A-WT_4,,hera,hera/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_hera.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20170918.A-WT_4,,StringTie,stringtie_onlyref/20170918.A-WT_4/20170918.A-WT_4.gtf,Rscripts/read_quant_stringtie.R,$(tx2geneext),yes))
$(eval $(call juncscalerule,20170918.A-WT_4,,SalmonMinimap2Nanopore,/home/Shared/data/seq/hussain_bath_nanopore_rnaseq/FGCZ/salmonminimap2/20171207_1645_p2557_4017_2_ALLREADS.pass/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),no))

$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,Salmon,salmon_stringtie_tx/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,SalmonBWA,salmonbwa_stringtie_tx/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,kallisto,kallisto_stringtie_tx/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_kallisto.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,RSEM,RSEM_stringtie_tx/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA.isoforms.results,Rscripts/read_quant_rsem.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,hera,hera_stringtie_tx/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_hera.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,StringTie,stringtie/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA_filtered.gtf,Rscripts/read_quant_stringtie.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes))

$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,Salmon,salmon_stringtie_tx/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,SalmonBWA,salmonbwa_stringtie_tx/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,kallisto,kallisto_stringtie_tx/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_kallisto.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,RSEM,RSEM_stringtie_tx/20170918.A-WT_4/20170918.A-WT_4.isoforms.results,Rscripts/read_quant_rsem.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,hera,hera_stringtie_tx/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_hera.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,StringTie,stringtie/20170918.A-WT_4/20170918.A-WT_4_filtered.gtf,Rscripts/read_quant_stringtie.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes))

## Combined coverages for all methods
define combcovrule
alpine/$(1)$(2)/alpine_combined_coverages.rds: STAR$(2)/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
alpine/$(1)$(2)/scaled_junction_coverage_Salmon.rds alpine/$(1)$(2)/scaled_junction_coverage_hera.rds \
alpine/$(1)$(2)/scaled_junction_coverage_RSEM.rds alpine/$(1)$(2)/scaled_junction_coverage_StringTie.rds \
alpine/$(1)$(2)/scaled_junction_coverage_SalmonBWA.rds alpine/$(1)$(2)/scaled_junction_coverage_kallisto.rds \
Rscripts/alpine_combine_scaled_coverages.R $(3) $(4)
	mkdir -p $$(@D)
	$(R) "--args junctioncovSTAR='STAR$(2)/$(1)/$(1)_SJ.out.tab' junctioncovSalmon='$$(word 2,$$^)' junctioncovSalmonBWA='$$(word 6,$$^)' junctioncovSalmonCDS='$(4)' junctioncovNanopore='$(3)' junctioncovhera='$$(word 3,$$^)' junctioncovkallisto='$$(word 7,$$^)' junctioncovRSEM='$$(word 4,$$^)' junctioncovStringTie='$$(word 5,$$^)' mmfracthreshold=$(mmfracthreshold) outrds='$$@'" Rscripts/alpine_combine_scaled_coverages.R Rout/alpine_combine_scaled_coverages_$(1)$(2).Rout
endef
$(eval $(call combcovrule,20151016.A-Cortex_RNA,,,alpine/20151016.A-Cortex_RNA/scaled_junction_coverage_SalmonCDS.rds))
$(eval $(call combcovrule,20170918.A-WT_4,,alpine/20170918.A-WT_4/scaled_junction_coverage_SalmonMinimap2Nanopore.rds,alpine/20170918.A-WT_4/scaled_junction_coverage_SalmonCDS.rds))
$(eval $(call combcovrule,20151016.A-Cortex_RNA,_stringtie_tx,,))
$(eval $(call combcovrule,20170918.A-WT_4,_stringtie_tx,,))

## Summarize gene expression from all methods
define combgexrule
alpine/$(1)$(2)/alpine_gene_expression.rds: alpine/$(1)$(2)/alpine_combined_coverages.rds \
Rscripts/combine_gene_expression_estimates.R
	$(R) "--args combcovrds='$$(word 1,$$^)' outrds='$$@'" Rscripts/combine_gene_expression_estimates.R Rout/combine_gene_expression_estimates_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call combgexrule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call combgexrule,$(notdir $(F)),_stringtie_tx)))

## Calculate gene score
define scorerule
alpine/$(1)$(2)/alpine_scores.rds: alpine/$(1)$(2)/alpine_combined_coverages.rds \
Rscripts/calculate_gene_scores.R
	$(R) "--args combcovrds='$$(word 1,$$^)' mmfracthreshold=$(mmfracthreshold) outrds='$$@'" Rscripts/calculate_gene_scores.R Rout/calculate_gene_scores_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call scorerule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call scorerule,$(notdir $(F)),_stringtie_tx)))



