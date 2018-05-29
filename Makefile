## Software paths
R := R_LIBS=/home/Shared/Rlib/devel-lib/ /usr/local/R/R-devel/bin/R CMD BATCH --no-restore --no-save
samtools := /usr/local/bin/samtools
bedtools := /usr/local/bin/bedtools

## List FASTQ files (without the _{R1,R2}.fastq.gz part)
fastqfilesreal := \
/home/Shared/data/seq/roche_pacbio_targeted_cdna/Illumina_RNA_seq/20151016.A-Cortex_RNA \
/home/Shared/data/seq/hussain_bath_nanopore_rnaseq/Illumina/FASTQ/20170918.A-WT_4
fastqfilessim := simulation/misannotated_utr/sim_misannotated_utr_1
fastqfiles := $(fastqfilesreal) $(fastqfilessim)

## Abundance quantification methods
quantmethods20151016.A-Cortex_RNA := Salmon SalmonSTAR kallisto RSEM StringTie hera SalmonCDS SalmonKeepDup
quantmethods20170918.A-WT_4 := Salmon SalmonSTAR kallisto RSEM StringTie hera SalmonCDS SalmonKeepDup SalmonMinimap2Nanopore WubMinimap2Nanopore
quantmethodssim_misannotated_utr_1 := Salmon SalmonSTAR kallisto RSEM StringTie hera SalmonCDS SalmonKeepDup
quantmethodsstringtie := Salmon SalmonSTAR kallisto RSEM StringTie hera

nthreads := 16

## Define the multimapping fraction threshold. Junctions with MM/(UM+MM)>mmfracthreshold will not be 
## included when calculating the "MM-aware" score
mmfracthreshold := 0.25

## Define the unique junction read threshold. Genes with less than this number of uniquely mapping 
## junction reads will be excluded from comparisons and evaluations
uniqjuncreadsthreshold := 25

## ==================================================================================== ##
##                                    Main rules                                        ##
## ==================================================================================== ##
.PHONY: all

all: prepref quant alpineprep scalecov \
preprefstringtie quantstringtie alpineprepstringtie scalecovstringtie plotsstringtie \
plots stats nanopore simulation

tmp: alpineprep quant scalecov plots

## Include makefiles. These need to be included after the definition of the "all" rule, 
## otherwise the default rule will be the first one from these makefiles
include makefiles/reference.mk
include makefiles/genes_to_plot.mk
include makefiles/Salmon.mk
include makefiles/hera.mk
include makefiles/SalmonSTAR.mk
include makefiles/StringTie.mk
include makefiles/kallisto.mk
include makefiles/RSEM.mk
include makefiles/STAR.mk
include makefiles/alpine.mk
include makefiles/plots.mk
include makefiles/stats.mk
include makefiles/simulate.mk
include makefiles/nanopore.mk

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
$(salmonkeepdupindex)/hash.bin \
$(kallistocdnancrnaindex) \
$(STARindextxome)/SA \
reference/RSEM/Homo_sapiens.GRCh38.rsem.cdna.ncrna/Homo_sapiens.GRCh38.rsem.cdna.ncrna.n2g.idx.fa \
$(STARindexnogtf)/chrNameLength.txt \
reference/hera/Homo_sapiens.GRCh38/index \
$(hisat2index).1.ht2 \
$(hisat2ss) \
$(gvizgenemodels)

## Align and quantify each sample
quant: $(foreach F,$(fastqfiles),salmon/cDNAncRNA/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),salmon/cDNAncRNAkeepdup/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),salmon/cds/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),kallisto/cDNAncRNA/$(notdir $(F))/abundance.tsv) \
$(foreach F,$(fastqfiles),salmonstartx/$(notdir $(F))/quant.sf) \
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
## CHESS annotation
########################################################################################################
## Prepare reference files and indexes
preprefchess: $(gtf_chess) $(txome_chess) $(tx2gene_chess) \
reference/salmon/chess2.0_assembly_fixed_sidx_v0.9.1/hash.bin \
reference/salmon/chess2.0_assembly_fixed_keepdup_sidx_v0.9.1/hash.bin \
reference/kallisto/chess2.0_assembly_fixed_kidx_v0.44.0 \
reference/STAR/chess2.0_assembly_fixed_STAR2.5.3a/SA \
reference/RSEM/chess2.0_assembly_fixed/chess2.0_assembly_fixed.n2g.idx.fa \
reference/hera/chess2.0_assembly_fixed/index \
reference/hisat2splicesites_chess2.0_assembly_fixed.txt

########################################################################################################
## Extended annotation (from StringTie)
########################################################################################################
## Prepare reference files and indexes
preprefstringtie: $(foreach F,$(fastqfiles),reference/$(notdir $(F))_stringtie_tx_tx2gene.rds) \
$(foreach F,$(fastqfiles),reference/$(notdir $(F))_stringtie_tx_tx2gene_withsymbol.rds) \
$(foreach F,$(fastqfiles),reference/salmon/$(notdir $(F))_stringtie_tx_sidx_v0.9.1/hash.bin) \
$(foreach F,$(fastqfiles),reference/kallisto/$(notdir $(F))_stringtie_tx_kidx_v0.44.0) \
$(foreach F,$(fastqfiles),reference/STAR/$(notdir $(F))_stringtie_tx_STAR2.5.3a/SA) \
$(foreach F,$(fastqfiles),reference/$(notdir $(F))_stringtie_tx_rsemgene2tx.txt) \
$(foreach F,$(fastqfiles),reference/RSEM/$(notdir $(F))_stringtie_tx/$(notdir $(F))_stringtie_tx.n2g.idx.fa) \
$(foreach F,$(fastqfiles),reference/hera/$(notdir $(F))_stringtie_tx/index) \
$(foreach F,$(fastqfiles),stringtie/$(notdir $(F))/$(notdir $(F))_filtered.gtf) \
$(foreach F,$(fastqfiles),stringtie/$(notdir $(F))/$(notdir $(F))_filtered_withgene.gtf) \
$(foreach F,$(fastqfiles),stringtie/$(notdir $(F))/$(notdir $(F))_stringtie_tx.fa)

## Align and quantify each sample
quantstringtie: $(foreach F,$(fastqfiles),salmon_stringtie_tx/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),kallisto_stringtie_tx/$(notdir $(F))/abundance.tsv) \
$(foreach F,$(fastqfiles),salmonstartx_stringtie_tx/$(notdir $(F))/quant.sf) \
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
## nanopore
########################################################################################################
nanopore: minimap2wub/$(nanopore_sample)/bam_count_reads.tsv \
minimap2genome/$(nanopore_sample)/$(nanopore_sample)_minimap2_genome_s.bam.bai \
minimap2genomebigwig/$(nanopore_sample)_minimap2_genome_s.bw \
minimap2salmon/$(nanopore_sample)/quant.sf

########################################################################################################
## Plots
########################################################################################################
plots: $(foreach F,$(fastqfiles),figures/observed_vs_predicted_junction_coverage/observed_vs_predicted_junction_coverage_$(notdir $(F)).rds) \
figures/predicted_coverage_pattern_comparison/predicted_coverage_pattern_comparison_20151016.A-Cortex_RNA_20170918.A-WT_4.rds \
$(foreach F,$(fastqfiles),figures/gene_scores/gene_scores_$(notdir $(F)).rds) \
$(foreach G,$(genes_to_plot),$(foreach F,$(fastqfiles),output_genewise/$(notdir $(F))/check/$(G).rds)) \
$(foreach F,$(fastqfiles),figures/correlation_with_inferential_variance/correlation_with_inferential_variance_$(notdir $(F)).rds) \
$(foreach F,$(fastqfiles),figures/correlation_between_methods/correlation_between_methods_$(notdir $(F)).rds) \
figures/correlation_between_nanopore_and_illumina_scores/correlation_between_nanopore_and_illumina_scores_20170918.A-WT_4.rds \
$(foreach F,$(fastqfiles),figures/correlation_with_true_abundances/correlation_with_true_abundances_$(notdir $(F)).rds) \
$(foreach F,$(fastqfiles),figures/association_exoncdscorrelation_score/association_exoncdscorrelation_score_$(notdir $(F)).rds) \
$(foreach F,$(fastqfilessim),figures/performance_simulated_data/performance_simulated_data_$(notdir $(F)).rds)

plotsstringtie: $(foreach F,$(fastqfiles),figures/observed_vs_predicted_junction_coverage/observed_vs_predicted_junction_coverage_$(notdir $(F))_stringtie_tx.rds) \
$(foreach F,$(fastqfiles),figures/gene_scores/gene_scores_$(notdir $(F))_stringtie_tx.rds)

########################################################################################################
## Stats
########################################################################################################
stats: $(foreach F,$(fastqfiles),stats/alpine_coverage_prediction_summary_$(notdir $(F)).txt) \
$(foreach F,$(fastqfiles),stats/alpine_coverage_prediction_summary_$(notdir $(F))_stringtie_tx.txt) \
$(foreach F,$(fastqfiles),stats/genes_with_high_score_$(notdir $(F)).txt) \
$(foreach F,$(fastqfiles),stats/genes_with_high_score_$(notdir $(F))_stringtie_tx.txt)

########################################################################################################
## Other
########################################################################################################
## List all packages
listpackages:
	$(R) Rscripts/list_packages.R Rout/list_packages.Rout



