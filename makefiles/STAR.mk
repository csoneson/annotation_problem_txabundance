STAR := /home/charlotte/software/STAR-2.5.3a/bin/Linux_x86_64/STAR
featurecounts := /home/charlotte/software/subread-1.6.0-Linux-x86_64/bin/featureCounts

STARindexnogtf := $(refdir)/genome/Homo_sapiens.GRCh38.dna.primary_assembly_STAR2.5.3a_nogtf

## STAR index without annotation gtf
$(STARindexnogtf)/SA: $(genome)
	mkdir -p $(@D)
	$(STAR) --runMode genomeGenerate --runThreadN $(nthreads) --genomeDir $(STARindexnogtf) \
	--genomeFastaFiles $(genome)

$(STARindexnogtf)/chrNameLength.txt: $(STARindexnogtf)/SA
	touch $(STARindexnogtf)/chrNameLength.txt

## ==================================================================================== ##
##                                     STAR                                             ##
## ==================================================================================== ##
## Align to the genome with STAR
define starrule
STAR$(3)/$(2)/$(2)_Aligned.sortedByCoord.out.bam: $(STARindexnogtf)/SA $(4) $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(STAR) --genomeDir $(STARindexnogtf) \
	--readFilesIn $(1)_R1.fastq.gz $(1)_R2.fastq.gz \
	--runThreadN $(nthreads) --outFileNamePrefix STAR$(3)/$(2)/$(2)_ \
	--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c \
	--sjdbGTFfile $(4) --sjdbOverhang $(5)
endef
$(foreach F,$(fastqfiles),$(eval $(call starrule,$(F),$(notdir $(F)),,$(gtf),150)))
$(foreach F,$(fastqfiles),$(eval $(call starrule,$(F),$(notdir $(F)),_stringtie_tx,stringtie/$(notdir $(F))/$(notdir $(F))_filtered.gtf,150)))

## Index STAR bam file
define starindexrule
STAR$(2)/$(1)/$(1)_Aligned.sortedByCoord.out.bam.bai: \
STAR$(2)/$(1)/$(1)_Aligned.sortedByCoord.out.bam
	$(samtools) index $$<
endef
$(foreach F,$(fastqfiles),$(eval $(call starindexrule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call starindexrule,$(notdir $(F)),_stringtie_tx)))

## Count reads mapping to exons and introns with featureCounts
define featurecountsrule
featureCounts$(3)/$(2)/$(2)_STAR_$(4).txt: \
STAR$(3)/$(2)/$(2)_Aligned.sortedByCoord.out.bam $(flatgtf$(4))
	mkdir -p $$(@D)
	$(featurecounts) -F GTF -t exon -g gene_id -O -s $(5) -p -T $(nthreads) -a $$(word 2,$$^) -o $$@ $$(word 1,$$^)
endef
$(foreach F,$(fastqfiles),$(eval $(call featurecountsrule,$(F),$(notdir $(F)),,exons,2)))
$(foreach F,$(fastqfiles),$(eval $(call featurecountsrule,$(F),$(notdir $(F)),,introns,2)))

## Convert BAM files to bigWig
define bwrule
STARbigwig/$(1)_Aligned.sortedByCoord.out.bw: $(STARindexnogtf)/chrNameLength.txt \
STAR/$(1)/$(1)_Aligned.sortedByCoord.out.bam
	mkdir -p $$(@D)	
	$(bedtools) genomecov -split -ibam $$(word 2,$$^) \
	-bg > STARbigwig/$(1)_Aligned.sortedByCoord.out.bedGraph
	
	bedGraphToBigWig STARbigwig/$(1)_Aligned.sortedByCoord.out.bedGraph \
	$$(word 1,$$^) $$@
	
	rm -f STARbigwig/$(1)_Aligned.sortedByCoord.out.bedGraph
endef
$(foreach F,$(fastqfiles),$(eval $(call bwrule,$(notdir $(F)))))

