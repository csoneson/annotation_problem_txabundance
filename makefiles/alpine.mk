## Fit bias model
define fitbiasrule
alpine/$(1)$(2)/alpine_fitbiasmodel.rds: $(3) STAR$(11)/$(1)/$(1)_Aligned.sortedByCoord.out.bam.bai \
Rscripts/alpine_fitbiasmodel.R
	mkdir -p $$(@D)
	$(R) "--args gtf='$(3)' bam='STAR$(11)/$(1)/$(1)_Aligned.sortedByCoord.out.bam' readlength=$(4) minsize=$(5) maxsize=$(6) organism='$(7)' genomeVersion='$(8)' version=$(9) outdir='$$(@D)' subsample=$(10)" Rscripts/alpine_fitbiasmodel.R Rout/alpine_fitbiasmodel_$(1)$(2).Rout
endef
$(eval $(call fitbiasrule,20151016.A-Cortex_RNA,,$(gtf),126,100,300,Homo_sapiens,GRCh38,90,TRUE,))
$(eval $(call fitbiasrule,20170918.A-WT_4,,$(gtf),151,140,450,Homo_sapiens,GRCh38,90,TRUE,))
#$(eval $(call fitbiasrule,SRR7056167,,$(gtf),XXX,XXX,XXX,Homo_sapiens,GRCh38,90,TRUE,))
$(eval $(call fitbiasrule,sim_misannotated_utr_1,,$(gtf),125,230,370,Homo_sapiens,GRCh38,90,FALSE,))
$(eval $(call fitbiasrule,20151016.A-Cortex_RNA,_stringtie_tx,stringtie/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA_filtered.gtf,126,100,300,Homo_sapiens,GRCh38,90,TRUE,_stringtie_tx))
$(eval $(call fitbiasrule,20170918.A-WT_4,_stringtie_tx,stringtie/20170918.A-WT_4/20170918.A-WT_4_filtered.gtf,151,140,450,Homo_sapiens,GRCh38,90,TRUE,_stringtie_tx))
#$(eval $(call fitbiasrule,SRR7056167,_stringtie_tx,stringtie/SRR7056167/SRR7056167_filtered.gtf,XXX,XXX,XXX,Homo_sapiens,GRCh38,90,TRUE,_stringtie_tx))
$(eval $(call fitbiasrule,sim_misannotated_utr_1,_stringtie_tx,stringtie/sim_misannotated_utr_1/sim_misannotated_utr_1_filtered.gtf,125,230,370,Homo_sapiens,GRCh38,90,TRUE,_stringtie_tx))
$(eval $(call fitbiasrule,20151016.A-Cortex_RNA,_chess,$(gtf_chess),126,100,300,Homo_sapiens,GRCh38,90,TRUE,_chess))
$(eval $(call fitbiasrule,20170918.A-WT_4,_chess,$(gtf_chess),151,140,450,Homo_sapiens,GRCh38,90,TRUE,_chess))
#$(eval $(call fitbiasrule,SRR7056167,_chess,$(gtf_chess),XXX,XXX,XXX,Homo_sapiens,GRCh38,90,TRUE,_chess))
$(eval $(call fitbiasrule,20151016.A-Cortex_RNA,_longUTR_added,$(gtf_longutr_added),126,100,300,Homo_sapiens,GRCh38,90,TRUE,))
$(eval $(call fitbiasrule,20170918.A-WT_4,_longUTR_added,$(gtf_longutr_added),151,140,450,Homo_sapiens,GRCh38,90,TRUE,))

## Predict transcript and junction coverage profiles for all transcripts that have at least one 
## junction and are longer than the fragment length
define predcovrule
alpine/$(1)$(2)/alpine_predicted_coverage.rds: alpine/$(1)$(2)/alpine_fitbiasmodel.rds \
STAR$(3)/$(1)/$(1)_Aligned.sortedByCoord.out.bam Rscripts/alpine_get_predicted_coverage.R
	mkdir -p $$(@D)
	$(R) "--args bam='STAR$(3)/$(1)/$(1)_Aligned.sortedByCoord.out.bam' biasmodels='alpine/$(1)$(2)/alpine_fitbiasmodel.rds' ncores=$(nthreads) outrds='$$@'" Rscripts/alpine_get_predicted_coverage.R Rout/alpine_get_predicted_coverage_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call predcovrule,$(notdir $(F)),,)))
$(foreach F,$(fastqfiles),$(eval $(call predcovrule,$(notdir $(F)),_stringtie_tx,_stringtie_tx)))
$(foreach F,$(fastqfilesreal),$(eval $(call predcovrule,$(notdir $(F)),_chess,_chess)))
$(foreach F,$(fastqfilesreal),$(eval $(call predcovrule,$(notdir $(F)),_longUTR_added,)))

## Scale junction coverage by transcript abundance estimates for each method
define juncscalerule
alpine/$(1)$(2)/scaled_junction_coverage_$(3).rds: alpine/$(1)$(2)/alpine_predicted_coverage.rds \
$(4) $(5) $(6) Rscripts/alpine_scale_junction_coverage.R
	mkdir -p $$(@D)
	$(R) "--args predcovrds='$$(word 1,$$^)' txquants='$(4)' quantreadscript='$(5)' tx2gene='$(6)' strandspec='$(7)' permutecounts=$(8) method='$(3)' outrds='$$@'" Rscripts/alpine_scale_junction_coverage.R Rout/alpine_scale_junction_coverage_$(1)$(2)_$(3).Rout
endef
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,Salmon,salmon/cDNAncRNA/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,SalmonSTAR,salmonstartx/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,SalmonKeepDup,salmon/cDNAncRNAkeepdup/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,SalmonCDS,salmon/cds/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,kallisto,kallisto/cDNAncRNA/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_kallisto.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,RSEM,RSEM/cDNAncRNA/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA.isoforms.results,Rscripts/read_quant_rsem.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,hera,hera/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_hera.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,,StringTie,stringtie_onlyref/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA.gtf,Rscripts/read_quant_stringtie.R,$(tx2geneext),yes,FALSE))

$(eval $(call juncscalerule,20170918.A-WT_4,,Salmon,salmon/cDNAncRNA/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,,SalmonSTAR,salmonstartx/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,,SalmonKeepDup,salmon/cDNAncRNAkeepdup/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,,SalmonCDS,salmon/cds/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,,kallisto,kallisto/cDNAncRNA/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_kallisto.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,,RSEM,RSEM/cDNAncRNA/20170918.A-WT_4/20170918.A-WT_4.isoforms.results,Rscripts/read_quant_rsem.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,,hera,hera/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_hera.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,,StringTie,stringtie_onlyref/20170918.A-WT_4/20170918.A-WT_4.gtf,Rscripts/read_quant_stringtie.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,,SalmonMinimap2Nanopore,minimap2salmon/20171207_1645_p2557_4017_2_ALLREADS.pass/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),no,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,,WubMinimap2Nanopore,minimap2wub/20171207_1645_p2557_4017_2_ALLREADS.pass/bam_count_reads.tsv,Rscripts/read_quant_wub.R,$(tx2geneext),no,FALSE))

$(eval $(call juncscalerule,20170918.A-WT_4,,SalmonPermuted,salmon/cDNAncRNA/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,TRUE))
$(eval $(call juncscalerule,20170918.A-WT_4,,SalmonSTARPermuted,salmonstartx/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,TRUE))
$(eval $(call juncscalerule,20170918.A-WT_4,,SalmonKeepDupPermuted,salmon/cDNAncRNAkeepdup/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,TRUE))
$(eval $(call juncscalerule,20170918.A-WT_4,,SalmonCDSPermuted,salmon/cds/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,TRUE))
$(eval $(call juncscalerule,20170918.A-WT_4,,kallistoPermuted,kallisto/cDNAncRNA/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_kallisto.R,$(tx2geneext),yes,TRUE))
$(eval $(call juncscalerule,20170918.A-WT_4,,RSEMPermuted,RSEM/cDNAncRNA/20170918.A-WT_4/20170918.A-WT_4.isoforms.results,Rscripts/read_quant_rsem.R,$(tx2geneext),yes,TRUE))
$(eval $(call juncscalerule,20170918.A-WT_4,,heraPermuted,hera/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_hera.R,$(tx2geneext),yes,TRUE))
$(eval $(call juncscalerule,20170918.A-WT_4,,StringTiePermuted,stringtie_onlyref/20170918.A-WT_4/20170918.A-WT_4.gtf,Rscripts/read_quant_stringtie.R,$(tx2geneext),yes,TRUE))

$(eval $(call juncscalerule,sim_misannotated_utr_1,,Salmon,salmon/cDNAncRNA/sim_misannotated_utr_1/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,,SalmonSTAR,salmonstartx/sim_misannotated_utr_1/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,,SalmonKeepDup,salmon/cDNAncRNAkeepdup/sim_misannotated_utr_1/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,,SalmonCDS,salmon/cds/sim_misannotated_utr_1/quant.sf,Rscripts/read_quant_salmon.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,,kallisto,kallisto/cDNAncRNA/sim_misannotated_utr_1/abundance.tsv,Rscripts/read_quant_kallisto.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,,RSEM,RSEM/cDNAncRNA/sim_misannotated_utr_1/sim_misannotated_utr_1.isoforms.results,Rscripts/read_quant_rsem.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,,hera,hera/sim_misannotated_utr_1/abundance.tsv,Rscripts/read_quant_hera.R,$(tx2geneext),yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,,StringTie,stringtie_onlyref/sim_misannotated_utr_1/sim_misannotated_utr_1.gtf,Rscripts/read_quant_stringtie.R,$(tx2geneext),yes,FALSE))

$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,Salmon,salmon_stringtie_tx/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,SalmonSTAR,salmonstartx_stringtie_tx/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,kallisto,kallisto_stringtie_tx/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_kallisto.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,RSEM,RSEM_stringtie_tx/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA.isoforms.results,Rscripts/read_quant_rsem.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,hera,hera_stringtie_tx/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_hera.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_stringtie_tx,StringTie,stringtie/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA_filtered.gtf,Rscripts/read_quant_stringtie.R,reference/20151016.A-Cortex_RNA_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))

$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,Salmon,salmon_stringtie_tx/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,SalmonSTAR,salmonstartx_stringtie_tx/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,kallisto,kallisto_stringtie_tx/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_kallisto.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,RSEM,RSEM_stringtie_tx/20170918.A-WT_4/20170918.A-WT_4.isoforms.results,Rscripts/read_quant_rsem.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,hera,hera_stringtie_tx/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_hera.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_stringtie_tx,StringTie,stringtie/20170918.A-WT_4/20170918.A-WT_4_filtered.gtf,Rscripts/read_quant_stringtie.R,reference/20170918.A-WT_4_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))

$(eval $(call juncscalerule,sim_misannotated_utr_1,_stringtie_tx,Salmon,salmon_stringtie_tx/sim_misannotated_utr_1/quant.sf,Rscripts/read_quant_salmon.R,reference/sim_misannotated_utr_1_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,_stringtie_tx,SalmonSTAR,salmonstartx_stringtie_tx/sim_misannotated_utr_1/quant.sf,Rscripts/read_quant_salmon.R,reference/sim_misannotated_utr_1_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,_stringtie_tx,kallisto,kallisto_stringtie_tx/sim_misannotated_utr_1/abundance.tsv,Rscripts/read_quant_kallisto.R,reference/sim_misannotated_utr_1_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,_stringtie_tx,RSEM,RSEM_stringtie_tx/sim_misannotated_utr_1/sim_misannotated_utr_1.isoforms.results,Rscripts/read_quant_rsem.R,reference/sim_misannotated_utr_1_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,_stringtie_tx,hera,hera_stringtie_tx/sim_misannotated_utr_1/abundance.tsv,Rscripts/read_quant_hera.R,reference/sim_misannotated_utr_1_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))
$(eval $(call juncscalerule,sim_misannotated_utr_1,_stringtie_tx,StringTie,stringtie/sim_misannotated_utr_1/sim_misannotated_utr_1_filtered.gtf,Rscripts/read_quant_stringtie.R,reference/sim_misannotated_utr_1_stringtie_tx_tx2gene_withsymbol.rds,yes,FALSE))

$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_chess,Salmon,salmon_chess/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_chess,SalmonSTAR,salmonstartx_chess/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_chess,SalmonKeepDup,salmon_chesskeepdup/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_chess,kallisto,kallisto_chess/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_kallisto.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_chess,RSEM,RSEM_chess/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA.isoforms.results,Rscripts/read_quant_rsem.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_chess,hera,hera_chess/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_hera.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_chess,StringTie,stringtie_chess/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA.gtf,Rscripts/read_quant_stringtie.R,$(tx2gene_chess_withsymbol),yes,FALSE))

$(eval $(call juncscalerule,20170918.A-WT_4,_chess,Salmon,salmon_chess/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_chess,SalmonSTAR,salmonstartx_chess/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_chess,SalmonKeepDup,salmon_chesskepdup/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_chess,kallisto,kallisto_chess/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_kallisto.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_chess,RSEM,RSEM_chess/20170918.A-WT_4/20170918.A-WT_4.isoforms.results,Rscripts/read_quant_rsem.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_chess,hera,hera_chess/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_hera.R,$(tx2gene_chess_withsymbol),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_chess,StringTie,stringtie_chess/20170918.A-WT_4/20170918.A-WT_4.gtf,Rscripts/read_quant_stringtie.R,$(tx2gene_chess_withsymbol),yes,FALSE))

$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_longUTR_added,Salmon,salmon_longUTR_added/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2gene_longutr_added),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_longUTR_added,SalmonKeepDup,salmon_longUTR_addedkeepdup/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,$(tx2gene_longutr_added),yes,FALSE))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,_longUTR_added,kallisto,kallisto_longUTR_added/20151016.A-Cortex_RNA/abundance.tsv, Rscripts/read_quant_kallisto.R,$(tx2gene_longutr_added),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_longUTR_added,Salmon,salmon_longUTR_added/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2gene_longutr_added),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_longUTR_added,SalmonKeepDup,salmon_longUTR_addedkeepdup/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,$(tx2gene_longutr_added),yes,FALSE))
$(eval $(call juncscalerule,20170918.A-WT_4,_longUTR_added,kallisto,kallisto_longUTR_added/20170918.A-WT_4/abundance.tsv, Rscripts/read_quant_kallisto.R,$(tx2gene_longutr_added),yes,FALSE))

## Combine coverages for all methods
define combcovrule
output/$(1)$(2)_combined_coverages.rds: STAR$(2)/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
alpine/$(1)$(2)/scaled_junction_coverage_Salmon.rds alpine/$(1)$(2)/scaled_junction_coverage_hera.rds \
alpine/$(1)$(2)/scaled_junction_coverage_RSEM.rds alpine/$(1)$(2)/scaled_junction_coverage_StringTie.rds \
alpine/$(1)$(2)/scaled_junction_coverage_SalmonSTAR.rds alpine/$(1)$(2)/scaled_junction_coverage_kallisto.rds \
Rscripts/combine_scaled_coverages.R $(3) $(4) $(5) $(6) $(7) $(8)
	mkdir -p $$(@D)
	$(R) "--args junctioncovSTAR='STAR$(2)/$(1)/$(1)_SJ.out.tab' junctioncovSalmon='$$(word 2,$$^)' junctioncovSalmonSTAR='$$(word 6,$$^)' junctioncovSalmonKeepDup='$(5)' junctioncovSalmonCDS='$(4)' junctioncovNanopore='$(3)' junctioncovhera='$$(word 3,$$^)' junctioncovkallisto='$$(word 7,$$^)' junctioncovRSEM='$$(word 4,$$^)' junctioncovStringTie='$$(word 5,$$^)' genecharacteristics='$(6)' exoncountstxt='$(7)' introncountstxt='$(8)' outrds='$$@'" Rscripts/combine_scaled_coverages.R Rout/combine_scaled_coverages_$(1)$(2).Rout
endef
$(eval $(call combcovrule,20151016.A-Cortex_RNA,,,alpine/20151016.A-Cortex_RNA/scaled_junction_coverage_SalmonCDS.rds,alpine/20151016.A-Cortex_RNA/scaled_junction_coverage_SalmonKeepDup.rds,output/gene_characteristics.rds,featureCounts/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA_STAR_exons.txt,featureCounts/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA_STAR_introns.txt))
$(eval $(call combcovrule,20170918.A-WT_4,,alpine/20170918.A-WT_4/scaled_junction_coverage_WubMinimap2Nanopore.rds,alpine/20170918.A-WT_4/scaled_junction_coverage_SalmonCDS.rds,alpine/20170918.A-WT_4/scaled_junction_coverage_SalmonKeepDup.rds,output/gene_characteristics.rds,featureCounts/20170918.A-WT_4/20170918.A-WT_4_STAR_exons.txt,featureCounts/20170918.A-WT_4/20170918.A-WT_4_STAR_introns.txt))
$(eval $(call combcovrule,sim_misannotated_utr_1,,,alpine/sim_misannotated_utr_1/scaled_junction_coverage_SalmonCDS.rds,alpine/sim_misannotated_utr_1/scaled_junction_coverage_SalmonKeepDup.rds,output/gene_characteristics.rds,featureCounts/sim_misannotated_utr_1/sim_misannotated_utr_1_STAR_exons.txt,featureCounts/sim_misannotated_utr_1/sim_misannotated_utr_1_STAR_introns.txt))
$(eval $(call combcovrule,20151016.A-Cortex_RNA,_stringtie_tx,,,,,,))
$(eval $(call combcovrule,20170918.A-WT_4,_stringtie_tx,,,,,,))
$(eval $(call combcovrule,sim_misannotated_utr_1,_stringtie_tx,,,,,,))

define combcovpermrule
output/$(1)$(2)_combined_coverages_permuted.rds: STAR$(2)/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
alpine/$(1)$(2)/scaled_junction_coverage_SalmonPermuted.rds alpine/$(1)$(2)/scaled_junction_coverage_heraPermuted.rds \
alpine/$(1)$(2)/scaled_junction_coverage_RSEMPermuted.rds alpine/$(1)$(2)/scaled_junction_coverage_StringTiePermuted.rds \
alpine/$(1)$(2)/scaled_junction_coverage_SalmonSTARPermuted.rds alpine/$(1)$(2)/scaled_junction_coverage_kallistoPermuted.rds \
Rscripts/combine_scaled_coverages.R $(3) $(4) $(5) $(6) $(7) $(8)
	mkdir -p $$(@D)
	$(R) "--args junctioncovSTAR='STAR$(2)/$(1)/$(1)_SJ.out.tab' junctioncovSalmon='$$(word 2,$$^)' junctioncovSalmonSTAR='$$(word 6,$$^)' junctioncovSalmonKeepDup='$(5)' junctioncovSalmonCDS='$(4)' junctioncovNanopore='$(3)' junctioncovhera='$$(word 3,$$^)' junctioncovkallisto='$$(word 7,$$^)' junctioncovRSEM='$$(word 4,$$^)' junctioncovStringTie='$$(word 5,$$^)' genecharacteristics='$(6)' exoncountstxt='$(7)' introncountstxt='$(8)' outrds='$$@'" Rscripts/combine_scaled_coverages.R Rout/combine_scaled_coverages_$(1)$(2)_permuted.Rout
endef
$(eval $(call combcovpermrule,20170918.A-WT_4,,,alpine/20170918.A-WT_4/scaled_junction_coverage_SalmonCDSPermuted.rds,alpine/20170918.A-WT_4/scaled_junction_coverage_SalmonKeepDupPermuted.rds,output/gene_characteristics.rds,featureCounts/20170918.A-WT_4/20170918.A-WT_4_STAR_exons.txt,featureCounts/20170918.A-WT_4/20170918.A-WT_4_STAR_introns.txt))

## Combine coverages for all methods, CHESS
define combcovchessrule
output/$(1)$(2)_combined_coverages.rds: STAR$(2)/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
alpine/$(1)$(2)/scaled_junction_coverage_Salmon.rds alpine/$(1)$(2)/scaled_junction_coverage_kallisto.rds \
output/gene_characteristics_chess.rds Rscripts/combine_scaled_coverages.R
	mkdir -p $$(@D)
	$(R) "--args junctioncovSTAR='STAR$(2)/$(1)/$(1)_SJ.out.tab' junctioncovSalmon='$$(word 2,$$^)' junctioncovSalmonSTAR='' junctioncovSalmonKeepDup='' junctioncovSalmonCDS='' junctioncovNanopore='' junctioncovhera='' junctioncovkallisto='$$(word 3,$$^)' junctioncovRSEM='' junctioncovStringTie='' genecharacteristics='$$(word 4,$$^)' exoncountstxt='' introncountstxt='' outrds='$$@'" Rscripts/combine_scaled_coverages.R Rout/combine_scaled_coverages_$(1)$(2).Rout
endef
$(eval $(call combcovchessrule,20151016.A-Cortex_RNA,_chess))
$(eval $(call combcovchessrule,20170918.A-WT_4,_chess))

## Combine coverages for all methods, longUTRadded
define combcovlongutraddedrule
output/$(1)$(2)_combined_coverages.rds: STAR$(3)/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
alpine/$(1)$(2)/scaled_junction_coverage_Salmon.rds alpine/$(1)$(2)/scaled_junction_coverage_kallisto.rds \
alpine/$(1)$(2)/scaled_junction_coverage_SalmonKeepDup.rds \
output/gene_characteristics.rds Rscripts/combine_scaled_coverages.R
	mkdir -p $$(@D)
	$(R) "--args junctioncovSTAR='STAR$(3)/$(1)/$(1)_SJ.out.tab' junctioncovSalmon='$$(word 2,$$^)' junctioncovSalmonSTAR='' junctioncovSalmonKeepDup='$$(word 4,$$^)' junctioncovSalmonCDS='' junctioncovNanopore='' junctioncovhera='' junctioncovkallisto='$$(word 3,$$^)' junctioncovRSEM='' junctioncovStringTie='' genecharacteristics='$$(word 5,$$^)' exoncountstxt='' introncountstxt='' outrds='$$@'" Rscripts/combine_scaled_coverages.R Rout/combine_scaled_coverages_$(1)$(2).Rout
endef
$(eval $(call combcovlongutraddedrule,20151016.A-Cortex_RNA,_longUTR_added,))
$(eval $(call combcovlongutraddedrule,20170918.A-WT_4,_longUTR_added,))

## Calculate gene scores and add to summary table
define scorerule
output/$(1)$(2)_combined_coverages$(3)_with_scores.rds: output/$(1)$(2)_combined_coverages$(3).rds \
Rscripts/calculate_gene_scores.R
	$(R) "--args combcovrds='$$(word 1,$$^)' mmfracthreshold=$(mmfracthreshold) outrds='$$@'" Rscripts/calculate_gene_scores.R Rout/calculate_gene_scores_$(1)$(2)$(3).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call scorerule,$(notdir $(F)),,)))
$(foreach F,$(fastqfiles),$(eval $(call scorerule,$(notdir $(F)),_stringtie_tx,)))
$(foreach F,$(fastqfilesreal),$(eval $(call scorerule,$(notdir $(F)),_chess,)))
$(foreach F,$(fastqfilesreal),$(eval $(call scorerule,$(notdir $(F)),_longUTR_added,)))
$(eval $(call scorerule,20170918.A-WT_4,,_permuted))



