STAR := /home/charlotte/software/STAR-2.5.3a/bin/Linux_x86_64/STAR
featurecounts := /home/charlotte/software/subread-1.6.0-Linux-x86_64/bin/featureCounts

STARindexnogtf := $(refdir)/genome/Homo_sapiens.GRCh38.dna.primary_assembly_STAR2.5.3a_nogtf

## STAR index without annotation gtf
$(STARindexnogtf)/SA: $(genome)
	mkdir -p $(@D)
	$(STAR) --runMode genomeGenerate --runThreadN $(nthreads) --genomeDir $(STARindexnogtf) \
	--genomeFastaFiles $(genome)

$(STARindexnogtf)/chrNameLength.txt: $(STARindexnogtf)/SA
	touch $@

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
$(foreach F,$(fastqfilesreal),$(eval $(call starrule,$(F),$(notdir $(F)),_chess,$(gtf_chess),150)))

## Index STAR bam file
define starindexrule
STAR$(2)/$(1)/$(1)_Aligned.sortedByCoord.out.bam.bai: \
STAR$(2)/$(1)/$(1)_Aligned.sortedByCoord.out.bam
	$(samtools) index $$<
endef
$(foreach F,$(fastqfiles),$(eval $(call starindexrule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call starindexrule,$(notdir $(F)),_stringtie_tx)))
$(foreach F,$(fastqfilesreal),$(eval $(call starindexrule,$(notdir $(F)),_chess)))

## Count reads mapping to exons and introns with featureCounts
define featurecountsrule
featureCounts$(3)/$(2)/$(2)_STAR_$(4).txt: \
STAR$(3)/$(2)/$(2)_Aligned.sortedByCoord.out.bam $(flatgtf$(4))
	mkdir -p $$(@D)
	$(featurecounts) -F GTF -t exon -g gene_id -O -s $(5) -p -T $(nthreads) -a $$(word 2,$$^) -o $$@ $$(word 1,$$^)
endef
$(foreach F,$(fastqfiles),$(eval $(call featurecountsrule,$(F),$(notdir $(F)),,exons,2)))
$(foreach F,$(fastqfiles),$(eval $(call featurecountsrule,$(F),$(notdir $(F)),,introns,2)))

## Count junction reads with featureCounts
define starjunccountrule
featureCounts$(3)/$(2)/$(2)_STAR.txt.jcounts: STAR$(3)/$(2)/$(2)_Aligned.sortedByCoord.out.bam
	mkdir -p $$(@D)
	$(featurecounts) -F GTF -t exon -g gene_id -f --minOverlap 1 --primary -O --splitOnly -s 2 -J -G $(genome) \
	-p -T $(nthreads) -a $(gtf) -o featureCounts$(3)/$(2)/$(2)_STAR.txt $$(word 1,$$^)
endef
$(foreach F,$(fastqfiles),$(eval $(call starjunccountrule,$(F),$(notdir $(F)),)))

## Convert BAM files to bigWig
define bwrule
STARbigwig$(2)/$(1)_Aligned.sortedByCoord.out.bw: $(STARindexnogtf)/chrNameLength.txt \
STAR$(2)/$(1)/$(1)_Aligned.sortedByCoord.out.bam
	mkdir -p $$(@D)	
	$(bedtools) genomecov -split -ibam $$(word 2,$$^) \
	-bg > STARbigwig$(2)/$(1)_Aligned.sortedByCoord.out.bedGraph
	
	bedGraphToBigWig STARbigwig$(2)/$(1)_Aligned.sortedByCoord.out.bedGraph \
	$$(word 1,$$^) $$@
	
	rm -f STARbigwig$(2)/$(1)_Aligned.sortedByCoord.out.bedGraph
endef
$(foreach F,$(fastqfiles),$(eval $(call bwrule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call bwrule,$(notdir $(F)),_chess)))
$(foreach F,$(fastqfiles),$(eval $(call bwrule,$(notdir $(F)),_stringtie_tx)))
