## Software paths
R := R_LIBS=/home/Shared/Rlib/release-3.5-lib/ /usr/local/R/R-3.4.0/bin/R CMD BATCH --no-restore --no-save
salmon := /home/charlotte/software/Salmon-0.8.2_linux_x86_64/bin/salmon
kallisto := /home/charlotte/software/kallisto_linux-v0.43.1/kallisto
RSEM := /home/charlotte/software/RSEM-1.3.0
STAR := /home/Shared_penticton/software/STAR/source/STAR
samtools := /usr/local/bin/samtools
bedtools := /usr/local/bin/bedtools
hisat2 := /home/charlotte/software/hisat2-2.1.0
stringtie := /home/charlotte/software/stringtie-1.3.3b.Linux_x86_64/stringtie
bwa := /home/charlotte/software/bwa/bwa
hera := /home/charlotte/software/hera/build/hera
hera_build := /home/charlotte/software/hera/build/hera_build

## Reference files
refdir := /home/Shared/data/annotation/Human/Ensembl_GRCh38.90
cdna := $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all.fa
cds := $(refdir)/cds/Homo_sapiens.GRCh38.cds.all.fa
ncrna := $(refdir)/ncRNA/Homo_sapiens.GRCh38.ncrna.fa
genome := $(refdir)/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf := $(refdir)/gtf/Homo_sapiens.GRCh38.90.gtf
txome := reference/Homo_sapiens.GRCh38.cdna.ncrna.fa

## Indexes
STARindex := $(refdir)/genome/Homo_sapiens.GRCh38.dna.primary_assembly_STAR_sjdbOverlap150
hisat2index := $(refdir)/genome/hisat2idx/Homo_sapiens.GRCh38.dna.primary_assembly
salmoncdnancrnaindex := $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.ncrna_sidx_v0.8.2
salmoncdsindex := $(refdir)/cds/Homo_sapiens.GRCh38.cds.all_sidx_v0.8.2
kallistocdnancrnaindex := $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.ncrna.all_kidx
rsemcdnancrnaindex := reference/RSEM/Homo_sapiens.GRCh38.rsem.cdna.ncrna
bwacdnancrnaindex := $(txome).sa
heraindex := reference/hera/Homo_sapiens.GRCh38

## Other annotation files
hisat2ss := reference/hisat2splicesites.txt
tx2gene := reference/Homo_sapiens.GRCh38.90_tx2gene.rds
rsemgene2tx := reference/Homo_sapiens.GRCh38.90_gene2tx.txt

## List FASTQ files (without the _{R1,R2}.fastq.gz part)
fastqfiles := \
/home/Shared/data/seq/roche_pacbio_targeted_cdna/Illumina_RNA_seq/20151016.A-Cortex_RNA \
/home/Shared/data/seq/hussain_bath_nanopore_rnaseq/Illumina/FASTQ/20170918.A-WT_3

## ==================================================================================== ##
##                                    Main rules                                        ##
## ==================================================================================== ##
.PHONY: all

all: prepref quant \
$(foreach F,$(fastqfiles),output/$(notdir $(F))_cdna_vs_cds.rds) \
$(foreach F,$(fastqfiles),alpine_check/$(notdir $(F))/genes_to_run.txt.rds)

prepref: $(txome) $(salmoncdnancrnaindex)/hash.bin $(salmoncdsindex)/hash.bin \
$(kallistocdnancrnaindex) $(rsemcdnancrnaindex).n2g.idx.fa $(STARindex)/SA $(tx2gene) \
$(bwacdnancrnaindex)

quant: $(foreach F,$(fastqfiles),salmon/cDNAncRNA/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),salmon/cds/$(notdir $(F))/quant.sf) \
$(foreach F,$(fastqfiles),STAR/$(notdir $(F))/$(notdir $(F))_Aligned.sortedByCoord.out.bam.bai) \
$(foreach F,$(fastqfiles),STARbigwig/$(notdir $(F))_Aligned.sortedByCoord.out.bw) \
$(foreach F,$(fastqfiles),kallisto/cDNAncRNA/$(notdir $(F))/abundance.tsv) \
$(foreach F,$(fastqfiles),RSEM/cDNAncRNA/$(notdir $(F))/$(notdir $(F)).isoforms.results) \
$(foreach F,$(fastqfiles),stringtie/$(notdir $(F))/$(notdir $(F)).gtf) \
$(foreach F,$(fastqfiles),stringtie_onlyref/$(notdir $(F))/$(notdir $(F)).gtf) \
$(foreach F,$(fastqfiles),salmonbwa/cDNAncRNA/$(notdir $(F))/quant.sf)

hera: $(heraindex)/index $(foreach F,$(fastqfiles),hera/$(notdir $(F))/abundance.tsv)

alpineprep: $(foreach F,$(fastqfiles),alpine/$(notdir $(F))/alpine_fitbiasmodel.rds) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/genes_to_run.txt) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/alpine_genemodels.rds)

## subset_genes_to_run.txt is a manually created file, which can be used to test a few genes
sumsub: $(foreach F,$(fastqfiles),alpine_check/$(notdir $(F))/subset_genes_to_run.txt.rds)

## ==================================================================================== ##
##                                   Reference files                                    ##
## ==================================================================================== ##
$(txome): $(cdna) $(ncrna)
	cat $(cdna) $(ncrna) > $(txome)

## Generate tx2gene
$(tx2gene): $(cdna) $(cds) $(ncrna) Rscripts/generate_tx2gene.R
	$(R) "--args cdna='$(word 1,$^)' cds='$(word 2,$^)' ncrna='$(word 3,$^)' outrds='$@'" Rscripts/generate_tx2gene.R Rout/generate_tx2gene.Rout 

$(rsemgene2tx): $(tx2gene) Rscripts/generate_rsemgene2tx.R
	$(R) "--args tx2gene='$(tx2gene)' rsemgene2tx='$@'" Rscripts/generate_rsemgene2tx.R Rout/generate_rsemgene2tx.Rout

## Build Salmon index for cDNA/ncRNA and CDS sequences
$(salmoncdnancrnaindex)/hash.bin: $(txome)
	$(salmon) index -t $< -k 25 -i $(@D) -p 10 --type quasi

$(salmoncdsindex)/hash.bin: $(cds)
	$(salmon) index -t $< -k 25 -i $(@D) -p 10 --type quasi

## Build kallisto index for cDNA/ncRNA sequences
$(kallistocdnancrnaindex): $(txome)
	$(kallisto) index -i $@ -k 25 $<

## Build RSEM index
$(rsemcdnancrnaindex).n2g.idx.fa: $(txome) $(rsemgene2tx) 
	mkdir -p $(@D)
	$(RSEM)/rsem-prepare-reference --transcript-to-gene-map $(rsemgene2tx) --bowtie --bowtie-path /usr/bin \
	$(txome) $(rsemcdnancrnaindex)

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

## BWA index
$(bwacdnancrnaindex): $(txome)
	$(bwa) index $(txome)

## hera index
$(heraindex)/index: $(genome) $(gtf)
	mkdir -p $(@D)
	$(hera_build) --fasta $(genome) --gtf $(gtf) --outdir $(@D)/

## ==================================================================================== ##
##                                      HERA                                            ##
## ==================================================================================== ##
## Run hera
define herarule
hera/$(notdir $(1))/abundance.tsv: $(heraindex)/index $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(hera) quant -i $(heraindex) -o $$(@D) -w 1 -t 10 $(1)_R1.fastq.gz $(1)_R2.fastq.gz
endef
$(foreach F,$(fastqfiles),$(eval $(call herarule,$(F))))

## ==================================================================================== ##
##                                  BWA + Salmon                                        ##
## ==================================================================================== ##
## Run BWA
define bwarule
$(3)/$(notdir $(1))/$(notdir $(1)).bam: $(2) $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(bwa) mem $(txome) $(1)_R1.fastq.gz $(1)_R2.fastq.gz | $(samtools) view -b -F 0x0800 -@ 10 - > $$@
endef
$(foreach F,$(fastqfiles),$(eval $(call bwarule,$(F),$(bwacdnancrnaindex),bwa/cDNAncRNA)))

## Run Salmon in alignment-based mode
define salmonbwarule
$(2)/$(notdir $(1))/quant.sf: bwa/cDNAncRNA/$(notdir $(1))/$(notdir $(1)).bam $(txome)
	mkdir -p $$(@D)
	$(salmon) quant -l A -a $$(word 1,$$^) -t $(txome) -p 10 --seqBias --gcBias -o $$(@D) 
endef
$(foreach F,$(fastqfiles),$(eval $(call salmonbwarule,$(F),salmonbwa/cDNAncRNA)))

## ==================================================================================== ##
##                                    Salmon                                            ##
## ==================================================================================== ##
## Run Salmon
define salmonrule
$(3)/$(notdir $(1))/quant.sf: $(2)/hash.bin $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(salmon) quant -i $$(word 1,$$(^D)) -l A -p 10 -1 $$(word 2,$$^) -2 $$(word 3,$$^) -o $$(@D) --seqBias --gcBias --posBias
endef
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(salmoncdnancrnaindex),salmon/cDNAncRNA)))
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(salmoncdsindex),salmon/cds)))

## Compare Salmon quants from cDNA/ncRNA and CDS
define salmoncomprule
output/$(notdir $(1))_cdna_vs_cds.rds: salmon/cDNAncRNA/$(notdir $(1))/quant.sf salmon/cds/$(notdir $(1))/quant.sf \
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
	--known-splicesite-infile $(hisat2ss) --summary-file $$@_summary.txt | \
	$(samtools) view -b -@ 10 - | $(samtools) sort -o $$@ -T $(notdir $(1))tmp -@ 10 - 
endef
$(foreach F,$(fastqfiles),$(eval $(call hisat2rule,$(F))))

## Run StringTie
define stringtierule
stringtie/$(notdir $(1))/$(notdir $(1)).gtf: HISAT2/$(notdir $(1))/$(notdir $(1)).bam
	mkdir -p $$(@D)
	$(stringtie) $$< -o $$@	-p 10 -G $(gtf) -A $$@.tab
endef
$(foreach F,$(fastqfiles),$(eval $(call stringtierule,$(F))))

## Run StringTie without assembly of new transcripts
define stringtierefrule
stringtie_onlyref/$(notdir $(1))/$(notdir $(1)).gtf: HISAT2/$(notdir $(1))/$(notdir $(1)).bam
	mkdir -p $$(@D)
	$(stringtie) $$< -o $$@	-p 10 -G $(gtf) -e -A $$@.tab
endef
$(foreach F,$(fastqfiles),$(eval $(call stringtierefrule,$(F))))

## ==================================================================================== ##
##                                     RSEM                                             ##
## ==================================================================================== ##
## Run RSEM
define rsemrule
$(3)/$(notdir $(1))/$(notdir $(1)).isoforms.results: $(2).n2g.idx.fa $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	bash -c '$(RSEM)/rsem-calculate-expression -p 10 --bowtie-path /usr/bin --paired-end <(gunzip -c $(1)_R1.fastq.gz) <(gunzip -c $(1)_R2.fastq.gz) $(2) $(3)/$(notdir $(1))/$(notdir $(1))'
endef
$(foreach F,$(fastqfiles),$(eval $(call rsemrule,$(F),$(rsemcdnancrnaindex),RSEM/cDNAncRNA)))

## ==================================================================================== ##
##                                   kallisto                                           ##
## ==================================================================================== ##
## Run kallisto
define kallistorule
$(3)/$(notdir $(1))/abundance.tsv: $(2) $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(kallisto) quant -i $(2) -o $$(@D) --bias -t 10 $(1)_R1.fastq.gz $(1)_R2.fastq.gz
endef
$(foreach F,$(fastqfiles),$(eval $(call kallistorule,$(F),$(kallistocdnancrnaindex),kallisto/cDNAncRNA)))

## ==================================================================================== ##
##                                     STAR                                             ##
## ==================================================================================== ##
## Align to the genome with STAR
define starrule
STAR/$(notdir $(1))/$(notdir $(1))_Aligned.sortedByCoord.out.bam: $(STARindex)/SA $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(STAR) --genomeDir $(STARindex) \
	--readFilesIn $(1)_R1.fastq.gz $(1)_R2.fastq.gz \
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
	mkdir -p $$(@D)	
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
alpine/$(1)/gene_characteristics.rds: $(gtf) salmon/cDNAncRNA/$(1)/quant.sf \
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
salmon/cDNAncRNA/$(1)/quant.sf kallisto/cDNAncRNA/$(1)/abundance.tsv RSEM/cDNAncRNA/$(1)/$(1).isoforms.results \
stringtie_onlyref/$(1)/$(1).gtf salmonbwa/cDNAncRNA/$(1)/quant.sf hera/$(1)/abundance.tsv \
Rscripts/alpine_prepare_for_comparison.R
	mkdir -p $$(@D)
	$(R) "--args gtf='$(gtf)' junctioncov='STAR/$(1)/$(1)_SJ.out.tab' quantsf='$$(word 3,$$^)' quantsfbwa='$$(word 7,$$^)' quantsfnanopore='$(2)' abundancetsv='$$(word 4,$$^)' heratsv='$$(word 8,$$^)' isoformsresults='$$(word 5,$$^)' stringtiegtf='$$(word 6,$$^)' outrds='$$@'" Rscripts/alpine_prepare_for_comparison.R Rout/alpine_prepare_for_comparison_$(1).Rout
endef
$(eval $(call alpinerefrule,20151016.A-Cortex_RNA,))
$(eval $(call alpinerefrule,20170918.A-WT_3,/home/Shared/data/seq/hussain_bath_nanopore_rnaseq/NSK007/salmonminimap2/SS2_wt_1/quant.sf))

## Predict coverage and compare to observed junction coverage
## "gene" can be either a gene ID or a text file with a list of genes to investigate
define alpinepredrule
alpine_check/$(1)/$(notdir $(2)).rds: alpine/$(1)/alpine_fitbiasmodel.rds STAR/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
alpine/$(1)/alpine_genemodels.rds Rscripts/alpine_compare_coverage.R \
STARbigwig/$(1)_Aligned.sortedByCoord.out.bw $(2)
	mkdir -p $$(@D)
	mkdir -p alpine_out/$(1)/plots
	mkdir -p alpine_out/$(1)/jcov
	mkdir -p alpine_out/$(1)/tpm
	mkdir -p alpine_out/$(1)/count	
	$(R) "--args gene='$(2)' bam='$$(word 2,$$^)' bigwig='$$(word 5,$$^)' ncores=$(3) genemodels='$$(word 3,$$^)' biasmodels='$$(word 1,$$^)' outdir='alpine_out/$(1)' checkdir='$$(@D)'" Rscripts/alpine_compare_coverage.R Rout/alpine_compare_coverage_$(1)_$(notdir $(2)).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call alpinepredrule,$(notdir $(F)),alpine/$(notdir $(F))/genes_to_run.txt,25)))
$(foreach F,$(fastqfiles),$(eval $(call alpinepredrule,$(notdir $(F)),alpine/$(notdir $(F))/subset_genes_to_run.txt,25)))



