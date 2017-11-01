R := R_LIBS=/home/Shared/Rlib/release-3.5-lib/ /usr/local/R/R-3.4.0/bin/R CMD BATCH --no-restore --no-save
salmon := /home/charlotte/software/Salmon-0.8.2_linux_x86_64/bin/salmon
refdir := /home/Shared/data/annotation/Human/Ensembl_GRCh38.90
cdna := $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all.fa
cds := $(refdir)/cds/Homo_sapiens.GRCh38.cds.all.fa
ncrna := $(refdir)/ncRNA/Homo_sapiens.GRCh38.ncrna.fa
STAR := /home/Shared_penticton/software/STAR/source/STAR
STARindex := $(refdir)/genome/Homo_sapiens.GRCh38.dna.primary_assembly_STAR_sjdbOverlap125
genome := $(refdir)/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf := $(refdir)/gtf/Homo_sapiens.GRCh38.90.gtf
samtools := /usr/local/bin/samtools
bedtools := /usr/local/bin/bedtools

fastqdir := /home/Shared/data/seq/roche_pacbio_targeted_cdna/Illumina_RNA_seq
fastqname := 20151016.A-Cortex_RNA

plotgene := ENSG00000198888 ENSG00000198712 ENSG00000162783 ENSG00000116120 ENSG00000174684 \
ENSG00000151067 ENSG00000100764 ENSG00000198728 ENSG00000185591 ENSG00000125354 \
ENSG00000158258 ENSG00000196814 ENSG00000153823 \
ENSG00000000003 ENSG00000000419 ENSG00000001497 ENSG00000004478 ENSG00000004779 \
ENSG00000005100 ENSG00000001630 ENSG00000280831 ENSG00000281709 ENSG00000284487

.PHONY: all

all: $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all_sidx_v0.8.2/hash.bin \
$(refdir)/cds/Homo_sapiens.GRCh38.cds.all_sidx_v0.8.2/hash.bin \
salmon/cDNA/$(fastqname)/quant.sf salmon/cds/$(fastqname)/quant.sf \
reference/Homo_sapiens.GRCh38.90_tx2gene.rds \
output/$(fastqname)_cdna_vs_cds.rds $(STARindex)/SA \
STARbigwig/$(fastqname)_Aligned.sortedByCoord.out.bw \
STAR/$(fastqname)/$(fastqname)_Aligned.sortedByCoord.out.bam.bai \
$(foreach G,$(plotgene),alpine_out/$(G).rds)

tmp: STAR/$(fastqname)/$(fastqname)_Aligned.sortedByCoord.out.bam

## ==================================================================================== ##
##                                   Reference files                                    ##
## ==================================================================================== ##
## Generate tx2gene
reference/Homo_sapiens.GRCh38.90_tx2gene.rds: $(cdna) $(cds) $(ncrna) Rscripts/generate_tx2gene.R
	$(R) "--args cdna='$(word 1,$^)' cds='$(word 2,$^)' ncrna='$(word 3,$^)' outrds='$@'" Rscripts/generate_tx2gene.R Rout/generate_tx2gene.Rout 

## Build Salmon index for cDNA and CDS sequences
$(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all_sidx_v0.8.2/hash.bin: $(cdna)
	$(salmon) index -t $< -k 25 -i $(@D) -p 10 --type quasi

$(refdir)/cds/Homo_sapiens.GRCh38.cds.all_sidx_v0.8.2/hash.bin: $(cds)
	$(salmon) index -t $< -k 25 -i $(@D) -p 10 --type quasi

## Build genome index
$(STARindex)/SA: $(genome) $(gtf)
	$(STAR) --runMode genomeGenerate --runThreadN 24 --genomeDir $(STARindex) \
	--genomeFastaFiles $(genome) --sjdbGTFfile $(gtf) --sjdbOverhang 125

## Get chromosome lengths
reference/chromosome_lengths.txt: STAR/$(fastqname)/$(fastqname)_Aligned.sortedByCoord.out.bam
	$(samtools) view -H $< | grep '@SQ' | cut -f2,3 | sed -e 's/SN://' | sed -e 's/LN://' > $@

## ==================================================================================== ##
##                                    Salmon                                            ##
## ==================================================================================== ##
## Run Salmon
salmon/cDNA/$(fastqname)/quant.sf: $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all_sidx_v0.8.2/hash.bin \
$(fastqdir)/$(fastqname)_R1.fastq.gz $(fastqdir)/$(fastqname)_R2.fastq.gz
	mkdir -p $(@D)
	$(salmon) quant -i $(word 1,$(^D)) -l A -p 10 -1 $(word 2,$^) -2 $(word 3,$^) -o $(@D) --seqBias --gcBias

salmon/cds/$(fastqname)/quant.sf: $(refdir)/cds/Homo_sapiens.GRCh38.cds.all_sidx_v0.8.2/hash.bin \
$(fastqdir)/$(fastqname)_R1.fastq.gz $(fastqdir)/$(fastqname)_R2.fastq.gz
	mkdir -p $(@D)
	$(salmon) quant -i $(word 1,$(^D)) -l A -p 10 -1 $(word 2,$^) -2 $(word 3,$^) -o $(@D) --seqBias --gcBias

## Compare Salmon quants from cDNA and CDS
output/$(fastqname)_cdna_vs_cds.rds: salmon/cDNA/$(fastqname)/quant.sf salmon/cds/$(fastqname)/quant.sf \
reference/Homo_sapiens.GRCh38.90_tx2gene.rds $(gtf) \
STARbigwig/$(fastqname)_Aligned.sortedByCoord.out.bw Rscripts/compare_cdna_and_cds_quants.R
	$(R) "--args cdnaquant='$(word 1,$^)' cdsquant='$(word 2,$^)' tx2gene='$(word 3,$^)' gtffile='$(word 4,$^)' bwfile='$(word 5,$^)' outrds='$@'" Rscripts/compare_cdna_and_cds_quants.R Rout/$(fastqname)_compare_cdna_and_cds_quants.Rout

## ==================================================================================== ##
##                                     STAR                                             ##
## ==================================================================================== ##
## Align to the genome with STAR
STAR/$(fastqname)/$(fastqname)_Aligned.sortedByCoord.out.bam: $(STARindex)/SA \
$(fastqdir)/$(fastqname)_R1.fastq.gz $(fastqdir)/$(fastqname)_R2.fastq.gz
	mkdir -p STAR/$(fastqname)
	$(STAR) --genomeDir $(STARindex) \
	--readFilesIn $(word 2,$^) $(word 3,$^) \
	--runThreadN 24 --outFileNamePrefix STAR/$(fastqname)/$(fastqname)_ \
	--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c

## Index STAR bam file
STAR/$(fastqname)/$(fastqname)_Aligned.sortedByCoord.out.bam.bai: \
STAR/$(fastqname)/$(fastqname)_Aligned.sortedByCoord.out.bam
	$(samtools) index $<

## Convert BAM files to bigWig
STARbigwig/$(fastqname)_Aligned.sortedByCoord.out.bw: reference/chromosome_lengths.txt \
STAR/$(fastqname)/$(fastqname)_Aligned.sortedByCoord.out.bam
	mkdir -p STARbigwig	
	$(bedtools) genomecov -split -ibam $(word 2,$^) \
	-bg > STARbigwig/$(fastqname)_Aligned.sortedByCoord.out.bedGraph
	
	bedGraphToBigWig STARbigwig/$(fastqname)_Aligned.sortedByCoord.out.bedGraph \
	$(word 1,$^) $@
	
	rm -f STARbigwig/$(fastqname)_Aligned.sortedByCoord.out.bedGraph

## ==================================================================================== ##
##                                   alpine                                             ##
## ==================================================================================== ##
## Fit bias model
alpine/alpine_fitbiasmodel.rds: $(gtf) STAR/$(fastqname)/$(fastqname)_Aligned.sortedByCoord.out.bam \
Rscripts/alpine_fitbiasmodel.R
	$(R) "--args gtf='$(gtf)' bam='$(word 2,$^)' outdir='alpine'" Rscripts/alpine_fitbiasmodel.R Rout/alpine_fitbiasmodel.Rout

## Predict coverage and compare to observed junction coverage
define alpinepredrule
alpine_out/$(1).rds: alpine/alpine_fitbiasmodel.rds STAR/$$(fastqname)/$$(fastqname)_Aligned.sortedByCoord.out.bam \
salmon/cDNA/$(fastqname)/quant.sf Rscripts/alpine_compare_coverage.R \
STARbigwig/$(fastqname)_Aligned.sortedByCoord.out.bw
	mkdir -p $$(@D)
	$(R) "--args gene='$(1)' bam='$$(word 2,$$^)' bw='$$(word 5,$$^)' gtf='$$(gtf)' quantsf='$$(word 3,$$^)' junctioncov='STAR/$$(fastqname)/$$(fastqname)_SJ.out.tab' biasmodels='$$(word 1,$$^)' outrds='$$@'" Rscripts/alpine_compare_coverage.R Rout/alpine_compare_coverage_$(1).Rout
endef
$(foreach G,$(plotgene),$(eval $(call alpinepredrule,$(G))))
