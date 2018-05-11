## Software paths
R := R_LIBS=/home/Shared/Rlib/devel-lib/ /usr/local/R/R-devel/bin/R CMD BATCH --no-restore --no-save
samtools := /usr/local/bin/samtools
bedtools := /usr/local/bin/bedtools

gvizgenemodels := reference/Homo_sapiens.GRCh38.90_gviz_genemodels.rds

## List FASTQ files (without the _{R1,R2}.fastq.gz part)
fastqfiles := \
/home/Shared/data/seq/roche_pacbio_targeted_cdna/Illumina_RNA_seq/20151016.A-Cortex_RNA \
/home/Shared/data/seq/hussain_bath_nanopore_rnaseq/Illumina/FASTQ/20170918.A-WT_4

## Abundance quantification methods
quantmethods20151016.A-Cortex_RNA := Salmon SalmonBWA kallisto RSEM StringTie hera SalmonCDS
quantmethods20170918.A-WT_4 := Salmon SalmonBWA kallisto RSEM StringTie hera SalmonCDS SalmonMinimap2Nanopore
quantmethodsstringtie := Salmon SalmonBWA kallisto RSEM StringTie hera

nthreads := 24

## Define the multimapping fraction threshold. Junctions with MM/(UM+MM)>mmfracthreshold will not be 
## included when calculating the "MM-aware" score
mmfracthreshold := 0.25

## ==================================================================================== ##
##                                    Main rules                                        ##
## ==================================================================================== ##
.PHONY: all

all: prepref quant alpineprep scalecov \
preprefstringtie quantstringtie alpineprepstringtie scalecovstringtie \
plots stats

## Include makefiles. These need to be included after the definition of the "all" rule, 
## otherwise the default rule will be the first one from these makefiles
include makefiles/reference.mk
include makefiles/Salmon.mk
include makefiles/hera.mk
include makefiles/SalmonBWA.mk
include makefiles/StringTie.mk
include makefiles/kallisto.mk
include makefiles/RSEM.mk
include makefiles/STAR.mk
include makefiles/alpine.mk
include makefiles/plots.mk
include makefiles/stats.mk
include makefiles/simulate.mk

########################################################################################################
## Original annotation
########################################################################################################
## Prepare reference files and indexes
prepref: $(txome) \
$(tx2gene) \
$(tx2geneext) \
$(flatgtfexons) \
$(salmoncdnancrnaindex)/hash.bin \
$(salmoncdsindex)/hash.bin \
$(kallistocdnancrnaindex) \
reference/bwa/Homo_sapiens.GRCh38.cdna.ncrna/Homo_sapiens.GRCh38.cdna.ncrna.fa.sa \
reference/RSEM/Homo_sapiens.GRCh38.rsem.cdna.ncrna/Homo_sapiens.GRCh38.rsem.cdna.ncrna.n2g.idx.fa \
$(STARindexnogtf)/chrNameLength.txt \
reference/hera/Homo_sapiens.GRCh38/index \
$(hisat2index).1.ht2 \
$(hisat2ss)

## Align and quantify each sample
quant: $(foreach F,$(fastqfiles),salmon/cDNAncRNA/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),salmon/cds/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),kallisto/cDNAncRNA/$(notdir $(F))/abundance.tsv) \
$(foreach F,$(fastqfiles),salmonbwa/cDNAncRNA/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),RSEM/cDNAncRNA/$(notdir $(F))/$(notdir $(F)).isoforms.results) \
$(foreach F,$(fastqfiles),STAR/$(notdir $(F))/$(notdir $(F))_Aligned.sortedByCoord.out.bam.bai) \
$(foreach F,$(fastqfiles),STARbigwig/$(notdir $(F))_Aligned.sortedByCoord.out.bw) \
$(foreach F,$(fastqfiles),$(foreach T,exons introns,featureCounts/$(notdir $(F))/$(notdir $(F))_STAR_$(T).txt)) \
$(foreach F,$(fastqfiles),hera/$(notdir $(F))/abundance.tsv) \
$(foreach F,$(fastqfiles),stringtie/$(notdir $(F))/$(notdir $(F)).gtf) \
$(foreach F,$(fastqfiles),stringtie_onlyref/$(notdir $(F))/$(notdir $(F)).gtf)

## Fit bias models and predict junction coverages
alpineprep: $(foreach F,$(fastqfiles),alpine/$(notdir $(F))/alpine_fitbiasmodel.rds) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/alpine_predicted_coverage.rds)

## Scale coverage and calculate scores
scalecov: $(foreach F,$(fastqfiles),$(foreach M,$(quantmethods$(F)),alpine/$(notdir $(F))/scaled_junction_coverage_$(M).rds)) \
$(foreach F,$(fastqfiles),output/$(notdir $(F))_combined_coverages.rds) \
$(foreach F,$(fastqfiles),output/$(notdir $(F))_combined_coverages_with_scores.rds)

########################################################################################################
## Extended annotation (from StringTie)
########################################################################################################
## Prepare reference files and indexes
preprefstringtie: $(foreach F,$(fastqfiles),reference/$(notdir $(F))_stringtie_tx_tx2gene.rds) \
$(foreach F,$(fastqfiles),reference/$(notdir $(F))_stringtie_tx_tx2gene_withsymbol.rds) \
$(foreach F,$(fastqfiles),reference/salmon/$(notdir $(F))_stringtie_tx_sidx_v0.9.1/hash.bin) \
$(foreach F,$(fastqfiles),reference/kallisto/$(notdir $(F))_stringtie_tx_kidx_v0.44.0) \
$(foreach F,$(fastqfiles),reference/bwa/$(notdir $(F))_stringtie_tx/$(notdir $(F))_stringtie_tx.fa.sa) \
$(foreach F,$(fastqfiles),reference/$(notdir $(F))_stringtie_tx_rsemgene2tx.txt) \
$(foreach F,$(fastqfiles),reference/RSEM/$(notdir $(F))_stringtie_tx/$(notdir $(F))_stringtie_tx.n2g.idx.fa) \
$(foreach F,$(fastqfiles),reference/hera/$(notdir $(F))_stringtie_tx/index) \
$(foreach F,$(fastqfiles),stringtie/$(notdir $(F))/$(notdir $(F))_filtered.gtf) \
$(foreach F,$(fastqfiles),stringtie/$(notdir $(F))/$(notdir $(F))_filtered_withgene.gtf) \
$(foreach F,$(fastqfiles),stringtie/$(notdir $(F))/$(notdir $(F))_stringtie_tx.fa)

## Align and quantify each sample
quantstringtie: $(foreach F,$(fastqfiles),salmon_stringtie_tx/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),kallisto_stringtie_tx/$(notdir $(F))/abundance.tsv) \
$(foreach F,$(fastqfiles),salmonbwa_stringtie_tx/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),RSEM_stringtie_tx/$(notdir $(F))/$(notdir $(F)).isoforms.results) \
$(foreach F,$(fastqfiles),STAR_stringtie_tx/$(notdir $(F))/$(notdir $(F))_Aligned.sortedByCoord.out.bam.bai) \
$(foreach F,$(fastqfiles),hera_stringtie_tx/$(notdir $(F))/abundance.tsv)

## Fit bias models and predict junction coverages
alpineprepstringtie: $(foreach F,$(fastqfiles),alpine/$(notdir $(F))_stringtie_tx/alpine_fitbiasmodel.rds) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))_stringtie_tx/alpine_predicted_coverage.rds)

## Scale coverage and calculate scores
scalecovstringtie: $(foreach F,$(fastqfiles),$(foreach M,$(quantmethodsstringtie),alpine/$(notdir $(F))_stringtie_tx/scaled_junction_coverage_$(M).rds)) \
$(foreach F,$(fastqfiles),output/$(notdir $(F))_stringtie_tx_combined_coverages.rds) \
$(foreach F,$(fastqfiles),output/$(notdir $(F))_stringtie_tx_combined_coverages_with_scores.rds)

########################################################################################################
## Simulated data
########################################################################################################
simulation: simulation/misannotated_utr/sim_misannotated_utr_1_R1.fastq.gz \
simulation/misannotated_utr/sim_misannotated_utr_1_R2.fastq.gz

########################################################################################################
## Plots
########################################################################################################
plots: $(foreach F,$(fastqfiles),figures/observed_vs_predicted_junction_coverage/observed_vs_predicted_junction_coverage_$(notdir $(F)).rds) \
$(foreach F,$(fastqfiles),figures/observed_vs_predicted_junction_coverage/observed_vs_predicted_junction_coverage_$(notdir $(F))_stringtie_tx.rds) \
figures/predicted_coverage_pattern_comparison/predicted_coverage_pattern_comparison_20151016.A-Cortex_RNA_20170918.A-WT_4.rds \
$(foreach F,$(fastqfiles),figures/gene_scores/gene_scores_$(notdir $(F)).rds) \
$(foreach F,$(fastqfiles),figures/gene_scores/gene_scores_$(notdir $(F))_stringtie_tx.rds)

########################################################################################################
## Stats
########################################################################################################
stats: $(foreach F,$(fastqfiles),stats/alpine_coverage_prediction_summary_$(notdir $(F)).txt) \
$(foreach F,$(fastqfiles),stats/alpine_coverage_prediction_summary_$(notdir $(F))_stringtie_tx.txt)

########################################################################################################
## Other
########################################################################################################
## List all packages
listpackages:
	$(R) Rscripts/list_packages.R Rout/list_packages.Rout

## Gene models for Gviz
$(gvizgenemodels): $(gtf) Rscripts/generate_genemodels.R Rscripts/plot_tracks.R
	$(R) "--args gtf='$(gtf)' outrds='$@'" Rscripts/generate_genemodels.R Rout/generate_genemodels.Rout

## Gene models for Gviz
## TODO: Fix generate_genemodels.R (exon_id column doesn't exist in the StringTie gtf, but it has exon_number
define gvizgmrule
reference/Gviz/$(notdir $(1))_stringtie_tx_gviz_genemodels.rds: stringtie/$(notdir $(1))/$(notdir $(1)).gtf \
Rscripts/generate_genemodels.R Rscripts/plot_tracks.R
	mkdir -p $$(@D)
	$$(R) "--args gtf='stringtie/$$(notdir $(1))/$$(notdir $(1)).gtf' outrds='$$@'" Rscripts/generate_genemodels.R Rout/generate_genemodels_$$(notdir $(1)).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call gvizgmrule,$(F))))


