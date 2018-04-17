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
quantmethods := Salmon SalmonBWA kallisto RSEM StringTie hera SalmonCDS
quantmethods2 := $(quantmethods) SalmonMinimap2Nanopore
quantmethods3 := Salmon SalmonBWA kallisto RSEM StringTie hera

nthreads := 24

## ==================================================================================== ##
##                                    Main rules                                        ##
## ==================================================================================== ##
.PHONY: all

all: prepref quant \
$(foreach F,$(fastqfiles),output/$(notdir $(F))_cdna_vs_cds.rds) \
$(foreach F,$(fastqfiles),alpine_check/$(notdir $(F))/genes_to_run.txt.rds) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/$(notdir $(F))_gene_scores.rds)

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

startmp: $(STARindexnogtf)/SA \
$(foreach F,$(fastqfiles),STARbigwig/$(notdir $(F))_Aligned.sortedByCoord.out.bw) \
$(foreach F,$(fastqfiles),STAR/$(notdir $(F))/$(notdir $(F))_Aligned.sortedByCoord.out.bam.bai) \
$(foreach K,exons introns,$(foreach F,$(fastqfiles),featureCounts/$(notdir $(F))/$(notdir $(F))_STAR_$(K).txt)) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/alpine_fitbiasmodel.rds) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/alpine_predicted_coverage.rds)

salmontmp: $(foreach F,$(fastqfiles),salmon/cDNAncRNA/$(notdir $(F))/quant.sf)

## Prepare reference files and indexes
prepref: $(txome) $(salmoncdnancrnaindex)/hash.bin $(salmoncdsindex)/hash.bin \
$(kallistocdnancrnaindex) $(rsemcdnancrnaindex).n2g.idx.fa $(STARindex)/SA $(tx2gene) \
$(bwacdnancrnaindex) reference/hera/Homo_sapiens.GRCh38/index $(hisat2index).1.ht2 $(hisat2ss) $(rsemgene2tx)

## Align and quantify each sample
quant: $(foreach F,$(fastqfiles),salmon/cDNAncRNA/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),salmon/cds/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),STAR/$(notdir $(F))/$(notdir $(F))_Aligned.sortedByCoord.out.bam.bai) \
$(foreach F,$(fastqfiles),STARbigwig/$(notdir $(F))_Aligned.sortedByCoord.out.bw) \
$(foreach F,$(fastqfiles),kallisto/cDNAncRNA/$(notdir $(F))/abundance.tsv) \
$(foreach F,$(fastqfiles),RSEM/cDNAncRNA/$(notdir $(F))/$(notdir $(F)).isoforms.results) \
$(foreach F,$(fastqfiles),stringtie/$(notdir $(F))/$(notdir $(F)).gtf) \
$(foreach F,$(fastqfiles),stringtie_onlyref/$(notdir $(F))/$(notdir $(F)).gtf) \
$(foreach F,$(fastqfiles),salmonbwa/cDNAncRNA/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),hera/$(notdir $(F))/abundance.tsv)

## Prepare files for alpine
alpineprep: $(foreach F,$(fastqfiles),alpine/$(notdir $(F))/alpine_fitbiasmodel.rds) \
$(foreach F,$(fastqfiles),gene_selection/$(notdir $(F))/genes_to_run.txt) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/alpine_predicted_coverage.rds)

## Scale coverage
scalecov: $(foreach M,$(quantmethods),$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/scaled_junction_coverage_$(M).rds))

## subset_genes_to_run.txt is a manually created file, which can be used to test a few genes
sumsub: $(foreach F,$(fastqfiles),alpine_check/$(notdir $(F))/subset_genes_to_run.txt.rds)

## Run methods on extended annotation from StringTie
stringtietx: $(foreach F,$(fastqfiles),reference/$(notdir $(F))_stringtie_tx_tx2gene.rds) \
$(foreach F,$(fastqfiles),reference/$(notdir $(F))_stringtie_tx_rsemgene2tx.txt) \
$(foreach F,$(fastqfiles),reference/RSEM_stringtie_tx/$(notdir $(F))_stringtie_tx/$(notdir $(F))_stringtie_tx.n2g.idx.fa) \
$(foreach F,$(fastqfiles),reference/bwa/$(notdir $(F))_stringtie_tx/$(notdir $(F))_stringtie_tx.fa.sa) \
$(foreach F,$(fastqfiles),salmon_stringtie_tx/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),kallisto_stringtie_tx/$(notdir $(F))/abundance.tsv) \
$(foreach F,$(fastqfiles),STAR_stringtie_tx/$(notdir $(F))/$(notdir $(F))_Aligned.sortedByCoord.out.bam.bai) \
$(foreach F,$(fastqfiles),RSEM_stringtie_tx/$(notdir $(F))/$(notdir $(F)).isoforms.results) \
$(foreach F,$(fastqfiles),salmonbwa_stringtie_tx/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),reference/hera/$(notdir $(F))_stringtie_tx/index) \
$(foreach F,$(fastqfiles),hera_stringtie_tx/$(notdir $(F))/abundance.tsv) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))_stringtie_tx/alpine_predicted_coverage.rds) \
$(foreach M,$(quantmethods3),$(foreach F,$(fastqfiles),alpine/$(notdir $(F))_stringtie_tx/scaled_junction_coverage_$(M).rds)) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))_stringtie_tx/alpine_combined_coverages.rds) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))_stringtie_tx/alpine_gene_expression.rds) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))_stringtie_tx/$(notdir $(F))_stringtie_tx_gene_scores.rds)

tmp3: $(foreach F,$(fastqfiles),RSEM_stringtie_tx/$(notdir $(F))/$(notdir $(F)).isoforms.results) \
$(foreach F,$(fastqfiles),salmonbwa_stringtie_tx/$(notdir $(F))/quant.sf)

tmp2: $(foreach F,$(fastqfiles),STAR_stringtie_tx/$(notdir $(F))/$(notdir $(F))_Aligned.sortedByCoord.out.bam.bai)

tmp: $(foreach F,$(fastqfiles),alpine/$(notdir $(F))_stringtie_tx/alpine_predicted_coverage.rds)

tmp4: $(foreach F,$(fastqfiles),STAR/$(notdir $(F))/$(notdir $(F))_Aligned.sortedByCoord.out.bam.bai) \
$(foreach F,$(fastqfiles),STAR_stringtie_tx/$(notdir $(F))/$(notdir $(F))_Aligned.sortedByCoord.out.bam.bai) \
$(foreach F,$(fastqfiles),salmonbwa_stringtie_tx/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),salmon_stringtie_tx/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),kallisto_stringtie_tx/$(notdir $(F))/abundance.tsv)

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


