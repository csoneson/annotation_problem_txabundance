## Software paths
R := R_LIBS=/home/Shared/Rlib/release-3.5-lib/ /usr/local/R/R-3.4.0/bin/R CMD BATCH --no-restore --no-save
salmon := /home/charlotte/software/Salmon-0.9.1_linux_x86_64/bin/salmon
kallisto := /home/charlotte/software/kallisto_linux-v0.44.0/kallisto
RSEM := /home/charlotte/software/RSEM-1.3.0
STAR := /home/Shared_penticton/software/STAR/source/STAR
samtools := /usr/local/bin/samtools
bedtools := /usr/local/bin/bedtools
hisat2 := /home/charlotte/software/hisat2-2.1.0
stringtie := /home/charlotte/software/stringtie-1.3.3b.Linux_x86_64/stringtie
bwa := /home/charlotte/software/bwa/bwa
hera := /home/charlotte/software/hera/build/hera
hera_build := /home/charlotte/software/hera/build/hera_build
strawberry := /home/charlotte/software/strawberry
gffread := /home/charlotte/software/gffread-0.9.12.Linux_x86_64/gffread

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
STARindexnogtf := $(refdir)/genome/Homo_sapiens.GRCh38.dna.primary_assembly_STAR_nogtf
hisat2index := $(refdir)/genome/hisat2idx/Homo_sapiens.GRCh38.dna.primary_assembly
salmoncdnancrnaindex := $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.ncrna_sidx_v0.9.1
salmoncdsindex := $(refdir)/cds/Homo_sapiens.GRCh38.cds.all_sidx_v0.9.1
kallistocdnancrnaindex := $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.ncrna.all_kidx_v0.44.0
rsemcdnancrnaindex := reference/RSEM/Homo_sapiens.GRCh38.rsem.cdna.ncrna
bwacdnancrnaindex := $(txome).sa
heraindex := reference/hera/Homo_sapiens.GRCh38

## Other annotation files
hisat2ss := reference/hisat2splicesites.txt
tx2gene := reference/Homo_sapiens.GRCh38.90_tx2gene.rds
rsemgene2tx := reference/Homo_sapiens.GRCh38.90_gene2tx.txt
tx2geneext := reference/Homo_sapiens.GRCh38.90_tx2gene_ext.rds
gvizgenemodels := reference/Homo_sapiens.GRCh38.90_gviz_genemodels.rds

## List FASTQ files (without the _{R1,R2}.fastq.gz part)
fastqfiles := \
/home/Shared/data/seq/roche_pacbio_targeted_cdna/Illumina_RNA_seq/20151016.A-Cortex_RNA \
/home/Shared/data/seq/hussain_bath_nanopore_rnaseq/Illumina/FASTQ/20170918.A-WT_4

## Abundance quantification methods
quantmethods := Salmon SalmonBWA kallisto RSEM StringTie hera SalmonCDS
quantmethods2 := $(quantmethods) SalmonMinimap2Nanopore

## ==================================================================================== ##
##                                    Main rules                                        ##
## ==================================================================================== ##
.PHONY: all

all: prepref quant \
$(foreach F,$(fastqfiles),output/$(notdir $(F))_cdna_vs_cds.rds) \
$(foreach F,$(fastqfiles),alpine_check/$(notdir $(F))/genes_to_run.txt.rds) \
$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/$(notdir $(F))_gene_scores.rds)

## Prepare reference files and indexes
prepref: $(txome) $(salmoncdnancrnaindex)/hash.bin $(salmoncdsindex)/hash.bin \
$(kallistocdnancrnaindex) $(rsemcdnancrnaindex).n2g.idx.fa $(STARindex)/SA $(tx2gene) \
$(bwacdnancrnaindex) $(heraindex)/index $(hisat2index).1.ht2 $(hisat2ss) $(rsemgene2tx)

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

scalecov: $(foreach M,$(quantmethods),$(foreach F,$(fastqfiles),alpine/$(notdir $(F))/scaled_junction_coverage_$(M).rds))

## subset_genes_to_run.txt is a manually created file, which can be used to test a few genes
sumsub: $(foreach F,$(fastqfiles),alpine_check/$(notdir $(F))/subset_genes_to_run.txt.rds)

tmp: $(foreach F,$(fastqfiles),STAR_stringtie_tx/$(notdir $(F))/$(notdir $(F))_Aligned.sortedByCoord.out.bam) \
$(foreach F,$(fastqfiles),kallisto_stringtie_tx/$(notdir $(F))/abundance.tsv)

## ==================================================================================== ##
##                                   Reference files                                    ##
## ==================================================================================== ##
$(txome): $(cdna) $(ncrna)
	cat $(cdna) $(ncrna) > $(txome)

## Generate tx2gene
$(tx2gene): $(cdna) $(cds) $(ncrna) Rscripts/generate_tx2gene.R
	$(R) "--args cdna='$(word 1,$^)' cds='$(word 2,$^)' ncrna='$(word 3,$^)' outrds='$@'" Rscripts/generate_tx2gene.R Rout/generate_tx2gene.Rout 

$(tx2geneext): $(tx2gene) $(gtf) Rscripts/extend_tx2gene.R
	$(R) "--args tx2gene='$(tx2gene)' gtf='$(gtf)' outrds='$@'" Rscripts/extend_tx2gene.R Rout/extend_tx2gene.Rout

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
	mkdir -p $(STARindex)
	$(STAR) --runMode genomeGenerate --runThreadN 24 --genomeDir $(STARindex) \
	--genomeFastaFiles $(genome) --sjdbGTFfile $(gtf) --sjdbOverhang 150

$(STARindex)/chrNameLength.txt: $(STARindex)/SA

## STAR index without annotation gtf
$(STARindexnogtf)/SA: $(genome)
	mkdir -p $(STARindexnogtf)
	$(STAR) --runMode genomeGenerate --runThreadN 10 --genomeDir $(STARindexnogtf) \
	--genomeFastaFiles $(genome)

$(STARindexnogtf)/chrNameLength.txt: $(STARindexnogtf)/SA

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

## Gene models for Gviz
$(gvizgenemodels): $(gtf) Rscripts/generate_genemodels.R Rscripts/plot_tracks.R
	$(R) "--args gtf='$(gtf)' outrds='$@'" Rscripts/generate_genemodels.R Rout/generate_genemodels.Rout

## ==================================================================================== ##
##                     Reference files for extended annotation                          ##
## ==================================================================================== ##
## We run StringTie below with the option of detecting unannotated transcripts. 
## The output is a gtf file that we also convert to a transcript fasta file. 
## We then generate indexes for the other methods based on this annotation.
## Note that this annotation is different between the different input samples

## Generate tx2gene
## TODO: FIX
#$(tx2gene): $(cdna) $(cds) $(ncrna) Rscripts/generate_tx2gene.R
#	$(R) "--args cdna='$(word 1,$^)' cds='$(word 2,$^)' ncrna='$(word 3,$^)' outrds='$@'" Rscripts/generate_tx2gene.R Rout/generate_tx2gene.Rout 

## TODO: FIX
#$(rsemgene2tx): $(tx2gene) Rscripts/generate_rsemgene2tx.R
#	$(R) "--args tx2gene='$(tx2gene)' rsemgene2tx='$@'" Rscripts/generate_rsemgene2tx.R Rout/generate_rsemgene2tx.Rout

## Build Salmon index for cDNA/ncRNA and CDS sequences
define salmonindexrule
reference/salmon/$(notdir $(1))_stringtie_tx_sidx_v0.9.1/hash.bin: \
stringtie/$(notdir $(1))/$(notdir $(1))_stringtie_tx.fa
	$(salmon) index -t $$< -k 25 -i $$(@D) -p 10 --type quasi
endef
$(foreach F,$(fastqfiles),$(eval $(call salmonindexrule,$(F))))

## Build kallisto index for cDNA/ncRNA sequences
define kallistoindexrule
reference/kallisto/$(notdir $(1))_stringtie_tx_kidx_v0.44.0: \
stringtie/$(notdir $(1))/$(notdir $(1))_stringtie_tx.fa
	mkdir -p $$(@D)
	$(kallisto) index -i $$@ -k 25 $$<
endef
$(foreach F,$(fastqfiles),$(eval $(call kallistoindexrule,$(F))))

## Build RSEM index
define rsemindexrule
reference/RSEM_extended/$(notdir $(1))_stringtie_tx/$(notdir $(1))_stringtie_tx.n2g.idx.fa: \
stringtie/$(notdir $(1))/$(notdir $(1))_stringtie_tx.fa \
## TODO: Add appropriate rsemgene2tx as dependency and in the recipe
	mkdir -p $$(@D)
	$(RSEM)/rsem-prepare-reference --transcript-to-gene-map XXXXXXX --bowtie --bowtie-path /usr/bin \
	stringtie/$$(notdir $(1))/$$(notdir $(1))_stringtie_tx.fa \
	reference/RSEM_extended/$$(notdir $(1))_stringtie_tx/$$(notdir $(1))_stringtie_tx
endef
$(foreach F,$(fastqfiles),$(eval $(call rsemindexrule,$(F))))

## BWA index
define bwaindexrule
reference/bwa/$(notdir $(1))_stringtie_tx/$(notdir $(1))_stringtie_tx.fa.sa: \
stringtie/$(notdir $(1))/$(notdir $(1))_stringtie_tx.fa
	mkdir -p $$(@D)
	scp $$< reference/bwa/$$(notdir $(1))_stringtie_tx/$$(notdir $(1))_stringtie_tx.fa
	$(bwa) index reference/bwa/$$(notdir $(1))_stringtie_tx/$$(notdir $(1))_stringtie_tx.fa
endef
$(foreach F,$(fastqfiles),$(eval $(call bwaindexrule,$(F))))

## hera index
## TODO: FIX
#define heraindexrule
#reference/hera/$(notdir $(1))_stringtie_tx/index: $(genome) stringtie/$(notdir $(1))/$(notdir $(1)).gtf
#	mkdir -p $$(@D)
#	$(hera_build) --fasta $$(genome) --gtf stringtie/$$(notdir $(1))/$$(notdir $(1)).gtf --outdir $$(@D)/
#endef
#$(foreach F,$(fastqfiles),$(eval $(call heraindexrule,$(F))))

## Gene models for Gviz
## TODO: Fix generate_genemodels.R (exon_id column doesn't exist in the StringTie gtf, but it has exon_number
define gvizgmrule
reference/Gviz/$(notdir $(1))_stringtie_tx_gviz_genemodels.rds: stringtie/$(notdir $(1))/$(notdir $(1)).gtf \
Rscripts/generate_genemodels.R Rscripts/plot_tracks.R
	mkdir -p $$(@D)
	$$(R) "--args gtf='stringtie/$$(notdir $(1))/$$(notdir $(1)).gtf' outrds='$$@'" Rscripts/generate_genemodels.R Rout/generate_genemodels_$$(notdir $(1)).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call gvizgmrule,$(F))))

## ==================================================================================== ##
##                                      HERA                                            ##
## ==================================================================================== ##
## Run hera
define herarule
hera/$(notdir $(1))/abundance.tsv: $(heraindex)/index# $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(hera) quant -i $(heraindex) -o $$(@D) -w 1 -t 10 $(1)_R1.fastq.gz $(1)_R2.fastq.gz
endef
$(foreach F,$(fastqfiles),$(eval $(call herarule,$(F))))

## ==================================================================================== ##
##                                  BWA + Salmon                                        ##
## ==================================================================================== ##
## Run BWA
define bwarule
$(3)/$(notdir $(1))/$(notdir $(1)).bam: $(2)# $(1)_R1.fastq.gz $(1)_R2.fastq.gz
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
$(3)/$(notdir $(1))/quant.sf: $(2)/hash.bin# $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(salmon) quant -i $$(word 1,$$(^D)) -l A -p 10 -1 $(1)_R1.fastq.gz -2 $(1)_R2.fastq.gz -o $$(@D) --seqBias --gcBias --posBias
endef
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(salmoncdnancrnaindex),salmon/cDNAncRNA)))
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(salmoncdsindex),salmon/cds)))
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),reference/salmon/$(notdir $(F))_stringtie_tx_sidx_v0.9.1,salmon_stringtie_tx)))

## Compare Salmon quants from cDNA/ncRNA and CDS
define salmoncomprule
output/$(notdir $(1))_cdna_vs_cds.rds: salmon/cDNAncRNA/$(notdir $(1))/quant.sf salmon/cds/$(notdir $(1))/quant.sf \
$(tx2gene) $(gtf) \
STARbigwig/$(notdir $(1))_Aligned.sortedByCoord.out.bw Rscripts/compare_cdna_and_cds_quants.R \
Rscripts/plot_tracks.R
	$(R) "--args cdnaquant='$$(word 1,$$^)' cdsquant='$$(word 2,$$^)' tx2gene='$$(word 3,$$^)' gtffile='$$(word 4,$$^)' bwfile='$$(word 5,$$^)' outrds='$$@'" Rscripts/compare_cdna_and_cds_quants.R Rout/compare_cdna_and_cds_quants_$(notdir $(1)).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call salmoncomprule,$(F))))

## ==================================================================================== ##
##                              HISAT2 + StringTie                                      ##
## ==================================================================================== ##
## Run HISAT2
define hisat2rule
HISAT2/$(notdir $(1))/$(notdir $(1)).bam: $(hisat2index).1.ht2 $(hisat2ss)# $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(hisat2)/hisat2 -p 10 -x $(hisat2index) --dta -1 $(1)_R1.fastq.gz -2 $(1)_R2.fastq.gz \
	--known-splicesite-infile $(hisat2ss) --summary-file $$@_summary.txt | \
	$(samtools) view -b -@ 10 - | $(samtools) sort -o $$@ -T $(notdir $(1))tmp -@ 10 - 
endef
$(foreach F,$(fastqfiles),$(eval $(call hisat2rule,$(F))))

## Run StringTie
define stringtierule
stringtie/$(notdir $(1))/$(notdir $(1)).gtf: HISAT2/$(notdir $(1))/$(notdir $(1)).bam $(gtf)
	mkdir -p $$(@D)
	$(stringtie) $$< -o $$@ -p 10 -G $(gtf) -A $$@.tab
endef
$(foreach F,$(fastqfiles),$(eval $(call stringtierule,$(F))))

## Run StringTie without assembly of new transcripts
define stringtierefrule
stringtie_onlyref/$(notdir $(1))/$(notdir $(1)).gtf: HISAT2/$(notdir $(1))/$(notdir $(1)).bam $(gtf)
	mkdir -p $$(@D)
	$(stringtie) $$< -o $$@ -p 10 -G $(gtf) -e -A $$@.tab
endef
$(foreach F,$(fastqfiles),$(eval $(call stringtierefrule,$(F))))

## Get transcript fasta from the StringTie gtf with assembled+reference features
define gffreadrule
stringtie/$(notdir $(1))/$(notdir $(1))_stringtie_tx.fa: stringtie/$(notdir $(1))/$(notdir $(1)).gtf \
$(genome)
	$(gffread) -w $$@ -g $(genome) $$<
endef
$(foreach F,$(fastqfiles),$(eval $(call gffreadrule,$(F))))

## ==================================================================================== ##
##                                     RSEM                                             ##
## ==================================================================================== ##
## Run RSEM
define rsemrule
$(3)/$(notdir $(1))/$(notdir $(1)).isoforms.results: $(2).n2g.idx.fa# $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	bash -c '$(RSEM)/rsem-calculate-expression -p 10 --bowtie-path /usr/bin --paired-end <(gunzip -c $(1)_R1.fastq.gz) <(gunzip -c $(1)_R2.fastq.gz) $(2) $(3)/$(notdir $(1))/$(notdir $(1))'
	rm -f $(3)/$(notdir $(1))/$(notdir $(1)).transcript.bam
endef
$(foreach F,$(fastqfiles),$(eval $(call rsemrule,$(F),$(rsemcdnancrnaindex),RSEM/cDNAncRNA)))

## ==================================================================================== ##
##                                   kallisto                                           ##
## ==================================================================================== ##
## Run kallisto
define kallistorule
$(3)/$(notdir $(1))/abundance.tsv: $(2)# $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(kallisto) quant -i $(2) -o $$(@D) --bias -t 10 $(1)_R1.fastq.gz $(1)_R2.fastq.gz
endef
$(foreach F,$(fastqfiles),$(eval $(call kallistorule,$(F),$(kallistocdnancrnaindex),kallisto/cDNAncRNA)))
$(foreach F,$(fastqfiles),$(eval $(call kallistorule,$(F),reference/kallisto/$(notdir $(F))_stringtie_tx_kidx_v0.44.0,kallisto_stringtie_tx)))

## ==================================================================================== ##
##                                     STAR                                             ##
## ==================================================================================== ##
## Align to the genome with STAR
define starrule
STAR/$(notdir $(1))/$(notdir $(1))_Aligned.sortedByCoord.out.bam: $(STARindex)/SA# $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(STAR) --genomeDir $(STARindex) \
	--readFilesIn $(1)_R1.fastq.gz $(1)_R2.fastq.gz \
	--runThreadN 10 --outFileNamePrefix STAR/$(notdir $(1))/$(notdir $(1))_ \
	--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c
endef
$(foreach F,$(fastqfiles),$(eval $(call starrule,$(F))))

## Align with StringTie-derived gtf files
define starrule2
STAR_stringtie_tx/$(notdir $(1))/$(notdir $(1))_Aligned.sortedByCoord.out.bam: $(STARindexnogtf)/SA $(2)
	mkdir -p $$(@D)
	$(STAR) --genomeDir $(STARindexnogtf) \
	--readFilesIn $(1)_R1.fastq.gz $(1)_R2.fastq.gz \
	--runThreadN 10 --outFileNamePrefix STAR_stringtie_tx/$(notdir $(1))/$(notdir $(1))_ \
	--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c \
	--sjdbGTFfile $(2) --sjdbOverhang $(3)
endef
$(foreach F,$(fastqfiles),$(eval $(call starrule2,$(F),stringtie/$(notdir $(F))/$(notdir $(F)).gtf,150)))

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
##                estimate and combine transcript/junction coverages                    ##
## ==================================================================================== ##
## Fit bias model
define fitbiasrule
alpine/$(1)/alpine_fitbiasmodel.rds: $(gtf) STAR/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
Rscripts/alpine_fitbiasmodel.R
	mkdir -p $$(@D)
	$(R) "--args gtf='$(gtf)' bam='$$(word 2,$$^)' readlength=$(2) minsize=$(3) maxsize=$(4) outdir='$$(@D)'" Rscripts/alpine_fitbiasmodel.R Rout/alpine_fitbiasmodel_$(1).Rout
endef
$(eval $(call fitbiasrule,20151016.A-Cortex_RNA,126,100,300))
$(eval $(call fitbiasrule,20170918.A-WT_4,151,140,450))

## Predict transcript and junction coverage profiles for all transcripts that have at least one 
## junction and are longer than the fragment length
define predcovrule
alpine/$(1)/alpine_predicted_coverage.rds: alpine/$(1)/alpine_fitbiasmodel.rds \
STAR/$(1)/$(1)_Aligned.sortedByCoord.out.bam Rscripts/alpine_get_predicted_coverage.R
	mkdir -p $$(@D)
	$(R) "--args bam='$$(word 2,$$^)' biasmodels='$$(word 1,$$^)' ncores=$(2) outrds='$$@'" Rscripts/alpine_get_predicted_coverage.R Rout/alpine_get_predicted_coverage_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call predcovrule,$(notdir $(F)),22)))

## Scale junction coverage by transcript abundance estimates for each method
define juncscalerule
alpine/$(1)/scaled_junction_coverage_$(2).rds: alpine/$(1)/alpine_predicted_coverage.rds \
$(3) $(4) $(tx2geneext) Rscripts/alpine_scale_junction_coverage.R
	mkdir -p $$(@D)
	$(R) "--args predcovrds='$$(word 1,$$^)' txquants='$(3)' quantreadscript='$(4)' tx2gene='$(tx2geneext)' strandspec='$(5)' method='$(2)' outrds='$$@'" Rscripts/alpine_scale_junction_coverage.R Rout/alpine_scale_junction_coverage_$(1)_$(2).Rout
endef
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,Salmon,salmon/cDNAncRNA/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,SalmonBWA,salmonbwa/cDNAncRNA/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,SalmonCDS,salmon/cds/20151016.A-Cortex_RNA/quant.sf,Rscripts/read_quant_salmon.R,yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,kallisto,kallisto/cDNAncRNA/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_kallisto.R,yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,RSEM,RSEM/cDNAncRNA/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA.isoforms.results,Rscripts/read_quant_rsem.R,yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,hera,hera/20151016.A-Cortex_RNA/abundance.tsv,Rscripts/read_quant_hera.R,yes))
$(eval $(call juncscalerule,20151016.A-Cortex_RNA,StringTie,stringtie_onlyref/20151016.A-Cortex_RNA/20151016.A-Cortex_RNA.gtf,Rscripts/read_quant_stringtie.R,yes))

$(eval $(call juncscalerule,20170918.A-WT_4,Salmon,salmon/cDNAncRNA/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,SalmonBWA,salmonbwa/cDNAncRNA/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,SalmonCDS,salmon/cds/20170918.A-WT_4/quant.sf,Rscripts/read_quant_salmon.R,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,kallisto,kallisto/cDNAncRNA/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_kallisto.R,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,RSEM,RSEM/cDNAncRNA/20170918.A-WT_4/20170918.A-WT_4.isoforms.results,Rscripts/read_quant_rsem.R,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,hera,hera/20170918.A-WT_4/abundance.tsv,Rscripts/read_quant_hera.R,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,StringTie,stringtie_onlyref/20170918.A-WT_4/20170918.A-WT_4.gtf,Rscripts/read_quant_stringtie.R,yes))
$(eval $(call juncscalerule,20170918.A-WT_4,SalmonMinimap2Nanopore,/home/Shared/data/seq/hussain_bath_nanopore_rnaseq/FGCZ/salmonminimap2/20171207_1645_p2557_4017_2_ALLREADS.pass/quant.sf,Rscripts/read_quant_salmon.R,no))

## Combined coverages for all methods
define combcovrule
alpine/$(1)/alpine_combined_coverages.rds: STAR/$(1)/$(1)_Aligned.sortedByCoord.out.bam \
alpine/$(1)/scaled_junction_coverage_Salmon.rds alpine/$(1)/scaled_junction_coverage_hera.rds \
alpine/$(1)/scaled_junction_coverage_RSEM.rds alpine/$(1)/scaled_junction_coverage_StringTie.rds \
alpine/$(1)/scaled_junction_coverage_SalmonBWA.rds alpine/$(1)/scaled_junction_coverage_kallisto.rds \
alpine/$(1)/scaled_junction_coverage_SalmonCDS.rds Rscripts/alpine_combine_scaled_coverages.R $(2)
	mkdir -p $$(@D)
	$(R) "--args junctioncovSTAR='STAR/$(1)/$(1)_SJ.out.tab' junctioncovSalmon='$$(word 2,$$^)' junctioncovSalmonBWA='$$(word 6,$$^)' junctioncovSalmonCDS='$$(word 8,$$^)' junctioncovNanopore='$(2)' junctioncovhera='$$(word 3,$$^)' junctioncovkallisto='$$(word 7,$$^)' junctioncovRSEM='$$(word 4,$$^)' junctioncovStringTie='$$(word 5,$$^)' outrds='$$@'" Rscripts/alpine_combine_scaled_coverages.R Rout/alpine_combine_scaled_coverages_$(1).Rout
endef
$(eval $(call combcovrule,20151016.A-Cortex_RNA,))
$(eval $(call combcovrule,20170918.A-WT_4,alpine/20170918.A-WT_4/scaled_junction_coverage_SalmonMinimap2Nanopore.rds))

## ==================================================================================== ##
##                            characterize genes                                        ##
## ==================================================================================== ##
## From annotation
output/characterize_genes.rds: $(gtf) $(txome) Rscripts/characterize_genes.R
	$(R) "--args gtf='$(gtf)' txome='$(txome)' outrds='$@'" Rscripts/characterize_genes.R Rout/characterize_genes.Rout

## Summarize gene expression from all methods
define combgexrule
alpine/$(1)/alpine_gene_expression.rds: alpine/$(1)/alpine_combined_coverages.rds \
Rscripts/combine_gene_expression_estimates.R
	$(R) "--args combcovrds='$$(word 1,$$^)' outrds='$$@'" Rscripts/combine_gene_expression_estimates.R Rout/combine_gene_expression_estimates_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call combgexrule,$(notdir $(F)))))

## ==================================================================================== ##
##                            plot gene scores                                          ##
## ==================================================================================== ##
define plotscorerule
alpine/$(1)/$(1)_gene_scores.rds: alpine/$(1)/alpine_gene_expression.rds \
output/characterize_genes.rds alpine/$(1)/alpine_combined_coverages.rds Rscripts/plot_score_distribution.R
	$(R) "--args covrds='$$(word 3,$$^)' gexrds='$$(word 1,$$^)' geneinfords='$$(word 2,$$^)' outrds='$$@'" Rscripts/plot_score_distribution.R Rout/plot_score_distribution_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call plotscorerule,$(notdir $(F)))))

## ==================================================================================== ##
##                            select genes and plot                                     ##
## ==================================================================================== ##
## Generate text file with genes to investigate further
define genelistrule
gene_selection/$(1)/genes_to_run.txt: alpine/$(1)/alpine_predicted_coverage.rds $(tx2geneext) Rscripts/list_genes_to_run.R
	$(R) "--args inrds='$$<' tx2gene='$$(word 2,$$^)' outtxt='$$@'" Rscripts/list_genes_to_run.R Rout/list_genes_to_run_$(1).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call genelistrule,$(notdir $(F)))))

## Predict coverage and compare to observed junction coverage
## "gene" can be either a gene ID or a text file with a list of genes to investigate
define alpinepredrule
alpine_check/$(1)/$(notdir $(2)).rds: $(gvizgenemodels) \
STARbigwig/$(1)_Aligned.sortedByCoord.out.bw alpine/$(1)/alpine_combined_coverages.rds \
$(2) Rscripts/alpine_compare_coverage.R Rscripts/plot_tracks.R
	mkdir -p $$(@D)
	mkdir -p alpine_out/$(1)/plots
	mkdir -p alpine_out/$(1)/jcov
	mkdir -p alpine_out/$(1)/tpm
	mkdir -p alpine_out/$(1)/count	
	$(R) "--args gene='$(2)' bigwig='$$(word 2,$$^)' ncores=$(3) genemodels='$$(word 1,$$^)' combcovrds='$$(word 3,$$^)' outdir='alpine_out/$(1)' libid='$(1)_' checkdir='$$(@D)'" Rscripts/alpine_compare_coverage.R Rout/alpine_compare_coverage_$(1)_$(notdir $(2)).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call alpinepredrule,$(notdir $(F)),gene_selection/$(notdir $(F))/genes_to_run.txt,25)))
$(foreach F,$(fastqfiles),$(eval $(call alpinepredrule,$(notdir $(F)),gene_selection/$(notdir $(F))/subset_genes_to_run.txt,25)))

