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

plotgenelist := genes_all.txt
# plotgene := \
# ENSG00000116120 ENSG00000174684 \
# ENSG00000151067 ENSG00000100764 ENSG00000198728 ENSG00000185591 ENSG00000125354 \
# ENSG00000158258 ENSG00000196814 ENSG00000153823 \
# ENSG00000000003 ENSG00000000419 ENSG00000001497 ENSG00000004478 ENSG00000004779 \
# ENSG00000005100 ENSG00000001630 ENSG00000182141 ENSG00000214013 ENSG00000136810 \
# ENSG00000186868 ENSG00000065183 \
# ENSG00000171161 ENSG00000073060 ENSG00000185658 ENSG00000121671 ENSG00000181163 \
# ENSG00000122733 ENSG00000120697 ENSG00000137992 ENSG00000188191 ENSG00000184271 \
# ENSG00000129245 ENSG00000172057 ENSG00000214022 ENSG00000134532 ENSG00000166340 \
# ENSG00000130717 ENSG00000119698 ENSG00000170100 ENSG00000105325 ENSG00000166224 \
# ENSG00000173258 ENSG00000124562 ENSG00000150768 ENSG00000198673 \
# ENSG00000144791 ENSG00000065923 ENSG00000099250 ENSG00000160963 \
# ENSG00000068097 ENSG00000122008 ENSG00000102078 ENSG00000144674 \
# ENSG00000185896 ENSG00000110046 ENSG00000171388 ENSG00000111052 \
# ENSG00000172977 ENSG00000214194 ENSG00000197162 ENSG00000171311 \
# ENSG00000183864 ENSG00000149932 ENSG00000144550 ENSG00000197183 \
# ENSG00000138166

.PHONY: all

all: $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all_sidx_v0.8.2/hash.bin \
$(refdir)/cds/Homo_sapiens.GRCh38.cds.all_sidx_v0.8.2/hash.bin \
salmon/cDNA/$(fastqname)/quant.sf salmon/cds/$(fastqname)/quant.sf \
reference/Homo_sapiens.GRCh38.90_tx2gene.rds \
output/$(fastqname)_cdna_vs_cds.rds $(STARindex)/SA \
STARbigwig/$(fastqname)_Aligned.sortedByCoord.out.bw \
STAR/$(fastqname)/$(fastqname)_Aligned.sortedByCoord.out.bam.bai \
$(foreach G,$(plotgenelist),alpine_check/$(G).rds)

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
## Summarize gene characteristics
alpine/gene_characteristics.rds: $(gtf) salmon/cDNA/$(fastqname)/quant.sf \
reference/Homo_sapiens.GRCh38.90_tx2gene.rds Rscripts/summarize_gene_characteristics.R
	$(R) "--args quantsf='$(word 2,$^)' gtf='$(word 1,$^)' tx2gene='$(word 3,$^)' outrds='$@'" Rscripts/summarize_gene_characteristics.R Rout/summarize_geene_characteristics.Rout

## Fit bias model
alpine/alpine_fitbiasmodel.rds: $(gtf) STAR/$(fastqname)/$(fastqname)_Aligned.sortedByCoord.out.bam \
Rscripts/alpine_fitbiasmodel.R
	$(R) "--args gtf='$(gtf)' bam='$(word 2,$^)' outdir='alpine'" Rscripts/alpine_fitbiasmodel.R Rout/alpine_fitbiasmodel.Rout

## Prepare reference files
alpine/alpine_genemodels.rds: $(gtf) STAR/$(fastqname)/$(fastqname)_Aligned.sortedByCoord.out.bam \
salmon/cDNA/$(fastqname)/quant.sf Rscripts/alpine_prepare_for_comparison.R
	$(R) "--args gtf='$(gtf)' junctioncov='STAR/$(fastqname)/$(fastqname)_SJ.out.tab' quantsf='$(word 3,$^)' outrds='$@'" Rscripts/alpine_prepare_for_comparison.R Rout/alpine_prepare_for_comparison.Rout

## Predict coverage and compare to observed junction coverage
## "gene" can be either a gene ID or a text file with a list of genes to investigate
define alpinepredrule
alpine_check/$(1).rds: alpine/alpine_fitbiasmodel.rds STAR/$$(fastqname)/$$(fastqname)_Aligned.sortedByCoord.out.bam \
alpine/alpine_genemodels.rds Rscripts/alpine_compare_coverage.R \
STARbigwig/$(fastqname)_Aligned.sortedByCoord.out.bw
	mkdir -p $$(@D)
	mkdir -p alpine_out
	$(R) "--args gene='$(1)' bam='$$(word 2,$$^)' bigwig='$$(word 5,$$^)' ncores=$(2) genemodels='$$(word 3,$$^)' biasmodels='$$(word 1,$$^)' outdir='alpine_out' checkdir='$$(@D)'" Rscripts/alpine_compare_coverage.R Rout/alpine_compare_coverage_$(1).Rout
endef
$(foreach G,$(plotgene),$(eval $(call alpinepredrule,$(G),1)))
$(foreach P,$(plotgenelist),$(eval $(call alpinepredrule,$(P),25)))