R := R_LIBS=/home/Shared/Rlib/release-3.5-lib/ /usr/local/R/R-3.4.0/bin/R CMD BATCH --no-restore --no-save
salmon := /home/charlotte/software/Salmon-0.8.2_linux_x86_64/bin/salmon
kallisto := /home/charlotte/software/kallisto_linux-v0.43.1/kallisto
RSEM := /home/charlotte/software/RSEM-1.3.0
refdir := /home/Shared/data/annotation/Human/Ensembl_GRCh38.90
cdna := $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all.fa
cds := $(refdir)/cds/Homo_sapiens.GRCh38.cds.all.fa
ncrna := $(refdir)/ncRNA/Homo_sapiens.GRCh38.ncrna.fa
STAR := /home/Shared_penticton/software/STAR/source/STAR
STARindex := $(refdir)/genome/Homo_sapiens.GRCh38.dna.primary_assembly_STAR_sjdbOverlap150
genome := $(refdir)/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf := $(refdir)/gtf/Homo_sapiens.GRCh38.90.gtf
samtools := /usr/local/bin/samtools
bedtools := /usr/local/bin/bedtools
tx2gene := reference/Homo_sapiens.GRCh38.90_tx2gene.rds
rsemgene2tx := reference/Homo_sapiens.GRCh38.90_gene2tx.txt
hisat2 := /home/charlotte/software/hisat2-2.1.0
hisat2index := $(refdir)/genome/hisat2idx/Homo_sapiens.GRCh38.dna.primary_assembly
hisat2ss := reference/hisat2splicesites.txt
stringtie := /home/charlotte/software/stringtie-1.3.3b.Linux_x86_64/stringtie

## List FASTQ files (without the _{R1,R2}.fastq.gz)
fastqfiles := \
/home/Shared/data/seq/roche_pacbio_targeted_cdna/Illumina_RNA_seq/20151016.A-Cortex_RNA \
/home/Shared/data/seq/hussain_bath_nanopore_rnaseq/Illumina/FASTQ/20170918.A-WT_3

.PHONY: all

all: quant \
$(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all_sidx_v0.8.2/hash.bin \
$(refdir)/cds/Homo_sapiens.GRCh38.cds.all_sidx_v0.8.2/hash.bin \
$(STARindex)/SA \
$(tx2gene) \
$(foreach F,$(fastqfiles),output/$(notdir $(F))_cdna_vs_cds.rds) \
$(foreach F,$(fastqfiles),alpine_check/$(notdir $(F))/genes_to_run.txt.rds)

quant: $(foreach F,$(fastqfiles),salmon/cDNA/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),salmon/cds/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),STAR/$(notdir $(F))/$(notdir $(F))_Aligned.sortedByCoord.out.bam.bai) \
$(foreach F,$(fastqfiles),STARbigwig/$(notdir $(F))_Aligned.sortedByCoord.out.bw) \
$(foreach F,$(fastqfiles),kallisto/cDNA/$(notdir $(F))/abundance.tsv) \
$(foreach F,$(fastqfiles),RSEM/cDNA/$(notdir $(F))/$(notdir $(F)).isoforms.results) \
$(foreach F,$(fastqfiles),HISAT2/$(notdir $(F))/$(notdir $(F)).bam)

alpineprep: $(foreach F,$(fastqfiles),alpine/$(notdir $(F))/alpine_fitbiasmodel.rds) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/genes_to_run.txt) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/alpine_genemodels.rds)

## ==================================================================================== ##
##                                   Reference files                                    ##
## ==================================================================================== ##
## Generate tx2gene
$(tx2gene): $(cdna) $(cds) $(ncrna) Rscripts/generate_tx2gene.R
	$(R) "--args cdna='$(word 1,$^)' cds='$(word 2,$^)' ncrna='$(word 3,$^)' outrds='$@'" Rscripts/generate_tx2gene.R Rout/generate_tx2gene.Rout 

$(rsemgene2tx): $(tx2gene) Rscripts/generate_rsemgene2tx.R
	$(R) "--args tx2gene='$(tx2gene)' rsemgene2tx='$@'" Rscripts/generate_rsemgene2tx.R Rout/generate_rsemgene2tx.Rout

## Build Salmon index for cDNA and CDS sequences
$(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all_sidx_v0.8.2/hash.bin: $(cdna)
	$(salmon) index -t $< -k 25 -i $(@D) -p 10 --type quasi

$(refdir)/cds/Homo_sapiens.GRCh38.cds.all_sidx_v0.8.2/hash.bin: $(cds)
	$(salmon) index -t $< -k 25 -i $(@D) -p 10 --type quasi

## Build kallisto index for cDNA sequences
$(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all_kidx: $(cdna)
	$(kallisto) index -i $@ -k 25 $<

## Build RSEM index
reference/RSEM/Homo_sapiens.GRCh38.rsem.n2g.idx.fa: $(cdna) $(rsemgene2tx) 
	mkdir -p $(@D)
	$(RSEM)/rsem-prepare-reference --transcript-to-gene-map $(rsemgene2tx) --bowtie --bowtie-path /usr/bin \
	$(cdna) reference/RSEM/Homo_sapiens.GRCh38.rsem

## Build genome index for STAR
$(STARindex)/SA: $(genome) $(gtf)
	$(STAR) --runMode genomeGenerate --runThreadN 24 --genomeDir $(STARindex) \
	--genomeFastaFiles $(genome) --sjdbGTFfile $(gtf) --sjdbOverhang 150

## Build genome index for HISAT2
$(hisat2index).1.ht2: $(genome)
	mkdir -p $(@D)
	$(hisat2)/hisat2-build -p 10 $(genome) $(hisat2index)

## Extract splice sites for HISAT2
$(hisat2ss): $(gtf)
	python $(hisat2)/hisat2_extract_splice_sites.py $(gtf) > $(hisat2ss)

## Extract list of junctions per transcript
reference/junctions_by_transcript.rds: $(gtf)
	$(R) "--args gtf='$(gtf)' outrds='$@'" Rscripts/get_junctions_per_transcript.R Rout/get_junctions_per_transcript.Rout

## ==================================================================================== ##
##                                    Salmon                                            ##
## ==================================================================================== ##
## Run Salmon
define salmonrule
$(3)/$(notdir $(1))/quant.sf: $(2) $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(salmon) quant -i $$(word 1,$$(^D)) -l A -p 10 -1 $$(word 2,$$^) -2 $$(word 3,$$^) -o $$(@D) --seqBias --gcBias
endef
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all_sidx_v0.8.2/hash.bin,salmon/cDNA)))
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(refdir)/cds/Homo_sapiens.GRCh38.cds.all_sidx_v0.8.2/hash.bin,salmon/cds)))

## Compare Salmon quants from cDNA and CDS
define salmoncomprule
output/$(notdir $(1))_cdna_vs_cds.rds: salmon/cDNA/$(notdir $(1))/quant.sf salmon/cds/$(notdir $(1))/quant.sf \
$(tx2gene) $(gtf) \
STARbigwig/$(notdir $(1))_Aligned.sortedByCoord.out.bw Rscripts/compare_cdna_and_cds_quants.R
	$(R) "--args cdnaquant='$$(word 1,$$^)' cdsquant='$$(word 2,$$^)' tx2gene='$$(word 3,$$^)' gtffile='$$(word 4,$$^)' bwfile='$$(word 5,$$^)' outrds='$$@'" Rscripts/compare_cdna_and_cds_quants.R Rout/$(notdir $(1))_compare_cdna_and_cds_quants.Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call salmoncomprule,$(F))))

## ==================================================================================== ##
##                              HISAT2 + StringTie                                      ##
## ==================================================================================== ##
## Run HISAT2
define hisat2rule
HISAT2/$(notdir $(1))/$(notdir $(1)).bam: $(hisat2index).1.ht2 $(1)_R1.fastq.gz $(1)_R2.fastq.gz $(hisat2ss)
	mkdir -p $$(@D)
	$(hisat2)/hisat2 -p 10 -x $(hisat2index) --dta -1 $(1)_R1.fastq.gz -2 $(1)_R2.fastq.gz \
	--known-splicesite-infile $(hisat2ss) | \
	$(samtools) view -b -@ 10 - | $(samtools) sort - $$@
endef
$(foreach F,$(fastqfiles),$(eval $(call hisat2rule,$(F))))

## Run StringTie
define stringtierule
stringtie/$(notdir $(1))/$(notdir $(1)).gtf: HISAT2/$(notdir $(1))/$(notdir $(1)).bam
	mkdir -p $$(@D)
	$(stringtie) $$< -o $$@	-p 10 -G $(gtf) -A $$@.tab
endef
$(foreach F,$(fastqfiles),$(eval $(call stringtierule,$(F))))

## ==================================================================================== ##
##                                     RSEM                                             ##
## ==================================================================================== ##
## Run RSEM
define rsemrule
$(3)/$(notdir $(1))/$(notdir $(1)).isoforms.results: $(2).n2g.idx.fa $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	bash -c '$(RSEM)/rsem-calculate-expression -p 10 --bowtie-path /usr/bin --paired-end <(gunzip -c $$(word 2,$$^)) <(gunzip -c $$(word 3,$$^)) $(2) $(3)/$(notdir $(1))/$(notdir $(1))'
endef
$(foreach F,$(fastqfiles),$(eval $(call rsemrule,$(F),reference/RSEM/Homo_sapiens.GRCh38.rsem,RSEM/cDNA)))

## ==================================================================================== ##
##                                   kallisto                                           ##
## ==================================================================================== ##
## Run kallisto
define kallistorule
$(3)/$(notdir $(1))/abundance.tsv: $(2) $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(kallisto) quant -i $$(word 1,$$^) -o $$(@D) --bias -t 10 $$(word 2,$$^) $$(word 3,$$^)
endef
$(foreach F,$(fastqfiles),$(eval $(call kallistorule,$(F),$(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all_kidx,kallisto/cDNA)))

## ==================================================================================== ##
##                                     STAR                                             ##
## ==================================================================================== ##
## Align to the genome with STAR
define starrule
STAR/$(notdir $(1))/$(notdir $(1))_Aligned.sortedByCoord.out.bam: $(STARindex)/SA $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p STAR/$(notdir $(1))
	$(STAR) --genomeDir $(STARindex) \
	--readFilesIn $$(word 2,$$^) $$(word 3,$$^) \
	--runThreadN 24 --outFileNamePrefix STAR/$(notdir $(1))/$(notdir $(1))_ \
	--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c
endef
$(foreach F,$(fastqfiles),$(eval $(call starrule,$(F))))

## Index STAR bam file
define starindexrule
STAR/$(notdir $(1))/$(notdir $(1))_Aligned.sortedByCoord.out.bam.bai: \
STAR/$(notdir $(1))/$(notdir $(1))_Aligned.sortedByCoord.out.bam
	$(samtools) index $$<
endef
$(foreach F,$(fastqfiles),$(eval $(call starindexrule,$(F))))

## Convert BAM files to bigWig
define bwrule
STARbigwig/$(notdir $(1))_Aligned.sortedByCoord.out.bw: $(STARindex)/chrNameLength.txt \
STAR/$(notdir $(1))/$(notdir $(1))_Aligned.sortedByCoord.out.bam
	mkdir -p STARbigwig	
	$(bedtools) genomecov -split -ibam $$(word 2,$$^) \
	-bg > STARbigwig/$(notdir $(1))_Aligned.sortedByCoord.out.bedGraph
	
	bedGraphToBigWig STARbigwig/$(notdir $(1))_Aligned.sortedByCoord.out.bedGraph \
	$$(word 1,$$^) $$@
	
	rm -f STARbigwig/$(notdir $(1))_Aligned.sortedByCoord.out.bedGraph
endef
$(foreach F,$(fastqfiles),$(eval $(call bwrule,$(F))))

## ==================================================================================== ##
##                                   alpine                                             ##
## ==================================================================================== ##
## Summarize gene characteristics
define genecharrule
alpine/$(1)/gene_characteristics.rds: $(gtf) salmon/cDNA/$(1)/quant.sf \
$(tx2gene) Rscripts/summarize_gene_characteristics.R
	mkdir -p $$(@D)
	$(R) "--args quantsf='$$(word 2,$$^)' gtf='$$(word 1,$$^)' tx2gene='$$(word 3,$$^)' outrds='$$@'" Rscripts/summarize_gene_characteristics.R Rout/summarize_gene_characteristics_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call genecharrule,$(notdir $(F)))))

## Generate text file with genes to investigate further
define genelistrule
alpine/$(1)/genes_to_run.txt: alpine/$(1)/gene_characteristics.rds
	$(R) "--args inrds='$$<' outtxt='$$@'" Rscripts/list_genes_to_run.R Rout/list_genes_to_run_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call genelistrule,$(notdir $(F)))))

## Fit bias model
define fitbiasrule
alpine/$(1)/alpine_fitbiasmodel.rds: $(gtf) STAR/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
Rscripts/alpine_fitbiasmodel.R
	mkdir -p $$(@D)
	$(R) "--args gtf='$(gtf)' bam='$$(word 2,$$^)' readlength=$(2) minsize=$(3) maxsize=$(4) outdir='$$(@D)'" Rscripts/alpine_fitbiasmodel.R Rout/alpine_fitbiasmodel_$(1).Rout
endef
$(eval $(call fitbiasrule,20151016.A-Cortex_RNA,126,100,300))
$(eval $(call fitbiasrule,20170918.A-WT_3,151,140,450))

## Prepare reference files
define alpinerefrule
alpine/$(1)/alpine_genemodels.rds: $(gtf) STAR/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
salmon/cDNA/$(1)/quant.sf Rscripts/alpine_prepare_for_comparison.R
	mkdir -p $$(@D)
	$(R) "--args gtf='$(gtf)' junctioncov='STAR/$(1)/$(1)_SJ.out.tab' quantsf='$$(word 3,$$^)' outrds='$$@'" Rscripts/alpine_prepare_for_comparison.R Rout/alpine_prepare_for_comparison_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call alpinerefrule,$(notdir $(F)))))

## Predict coverage and compare to observed junction coverage
## "gene" can be either a gene ID or a text file with a list of genes to investigate
define alpinepredrule
alpine_check/$(1)/$(notdir $(2)).rds: alpine/$(1)/alpine_fitbiasmodel.rds STAR/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
alpine/$(1)/alpine_genemodels.rds Rscripts/alpine_compare_coverage.R \
STARbigwig/$(1)_Aligned.sortedByCoord.out.bw $(2)
	mkdir -p $$(@D)
	mkdir -p alpine_out/$(1)
	$(R) "--args gene='$(2)' bam='$$(word 2,$$^)' bigwig='$$(word 5,$$^)' ncores=$(3) genemodels='$$(word 3,$$^)' biasmodels='$$(word 1,$$^)' outdir='alpine_out/$(1)' checkdir='$$(@D)'" Rscripts/alpine_compare_coverage.R Rout/alpine_compare_coverage_$(1)_$(notdir $(2)).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call alpinepredrule,$(notdir $(F)),alpine/$(notdir $(F))/genes_to_run.txt,25)))



