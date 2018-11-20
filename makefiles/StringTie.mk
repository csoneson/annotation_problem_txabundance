hisat2 := /home/charlotte/software/hisat2-2.1.0
stringtie := /home/charlotte/software/stringtie-1.3.3b.Linux_x86_64/stringtie
gffread := /home/charlotte/software/gffread-0.9.12.Linux_x86_64/gffread
featurecounts := /home/charlotte/software/subread-1.6.0-Linux-x86_64/bin/featureCounts

hisat2index := $(refdir)/genome/hisat2idx/Homo_sapiens.GRCh38.dna.primary_assembly
hisat2ss := reference/hisat2splicesites.txt

## Build genome index for HISAT2
$(hisat2index).1.ht2: $(genome)
	mkdir -p $(@D)
	$(hisat2)/hisat2-build -p $(nthreads) $(genome) $(hisat2index)

## Extract splice sites for HISAT2
$(hisat2ss): $(gtf)
	python $(hisat2)/hisat2_extract_splice_sites.py $(gtf) > $(hisat2ss)

reference/hisat2splicesites_chess2.0_assembly_fixed.txt: $(gtf_chess)
	python $(hisat2)/hisat2_extract_splice_sites.py $(gtf_chess) > $@

## ==================================================================================== ##
##                              HISAT2 + StringTie                                      ##
## ==================================================================================== ##
## Run HISAT2
define hisat2rule
HISAT2$(3)/$(2)/$(2).bam: $(hisat2index).1.ht2 $(4) $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(hisat2)/hisat2 -p $(nthreads) -x $(hisat2index) --dta -1 $(1)_R1.fastq.gz -2 $(1)_R2.fastq.gz \
	--known-splicesite-infile $(4) --summary-file $$@_summary.txt | \
	$(samtools) view -b -@ $(nthreads) - | $(samtools) sort -o $$@ -T $(2)tmp -@ $(nthreads) - 
endef
$(foreach F,$(fastqfiles),$(eval $(call hisat2rule,$(F),$(notdir $(F)),,$(hisat2ss))))
$(foreach F,$(fastqfilesreal),$(eval $(call hisat2rule,$(F),$(notdir $(F)),_chess,reference/hisat2splicesites_chess2.0_assembly_fixed.txt)))

## Get junction counts from HiSAT2 alignment
#define hisatjunccountrule
#HISAT2$(3)/$(2)/$(2).jcounts: HISAT2$(3)/$(2)/$(2).bam
#	$(featurecounts) 


## Run StringTie without assembly of new transcripts
define stringtierefrule
stringtie_onlyref$(3)/$(1)/$(1).gtf: HISAT2$(3)/$(1)/$(1).bam $(4)
	mkdir -p $$(@D)
	$(stringtie) $$< -o $$@ -p $(nthreads) -G $(4) -e -A $$@.tab $(2)
endef
$(foreach F,$(fastqfiles),$(eval $(call stringtierefrule,$(notdir $(F)),--rf,,$(gtf))))
$(foreach F,$(fastqfilesreal),$(eval $(call stringtierefrule,$(notdir $(F)),--rf,_chess,$(gtf_chess))))

## Run StringTie with assembly of new transcripts
define stringtierule
stringtie/$(1)/$(1).gtf: HISAT2/$(1)/$(1).bam $(gtf)
	mkdir -p $$(@D)
	$(stringtie) $$< -o $$@ -p $(nthreads) -G $(gtf) -A $$@.tab $(2)
endef
$(foreach F,$(fastqfiles),$(eval $(call stringtierule,$(notdir $(F)),--rf)))

## Filter StringTie output gtf. For assembled single-exon transcripts, StringTie can not 
## derive the strand (since there are no reads spanning junctions). This will prevent the 
## gtf from being read with e.g. GenomicFeatures::makeTxDbFromGff. Thus, here we filter 
## out such transcripts. We also change the transcript_id to the reference_id, and the 
## gene_id to ref_gene_id, whenever the latter are available (annotated genes/transcripts)
define stringtiefilterrule
stringtie/$(1)/$(1)_filtered.gtf: stringtie/$(1)/$(1).gtf \
Rscripts/filter_stringtie_gtf.R
	$(R) "--args ingtf='$$<' outgtf='$$@'" Rscripts/filter_stringtie_gtf.R Rout/filter_stringtie_gtf_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call stringtiefilterrule,$(notdir $(F)))))

## Hera requires that there are "gene" entries as well as "exon" and "transcript" in the gtf file
define stringtiegenerule
stringtie/$(1)/$(1)_filtered_withgene.gtf: stringtie/$(1)/$(1)_filtered.gtf Rscripts/add_gene_to_stringtie_gtf.R
	$(R) "--args ingtf='$$<' outgtf='$$@'" Rscripts/add_gene_to_stringtie_gtf.R Rout/add_gene_to_stringtie_gtf_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call stringtiegenerule,$(notdir $(F)))))

## Get transcript fasta from the StringTie gtf with assembled+reference features
define gffreadrule
stringtie/$(1)/$(1)_stringtie_tx.fa: stringtie/$(1)/$(1)_filtered.gtf \
$(genome)
	$(gffread) -w $$@ -g $(genome) $$<
endef
$(foreach F,$(fastqfiles),$(eval $(call gffreadrule,$(notdir $(F)))))

## Generate gene models for Gviz
define stringtiegvizrule
stringtie/$(1)/$(1)_gviz_genemodels.rds: stringtie/$(1)/$(1)_filtered.gtf Rscripts/generate_genemodels.R Rscripts/helper_plot_tracks.R
	mkdir -p $$(@D)
	$(R) "--args gtf='$$(word 1,$$^)' outrds='$$@'" Rscripts/generate_genemodels.R Rout/generate_genemodels_$(1)_stringtie_tx.Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call stringtiegvizrule,$(notdir $(F)))))

