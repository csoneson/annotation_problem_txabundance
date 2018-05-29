STAR := /home/charlotte/software/STAR-2.5.3a/bin/Linux_x86_64/STAR
salmon := /home/charlotte/software/Salmon-0.9.1_linux_x86_64/bin/salmon

STARindextxome := reference/STAR/Homo_sapiens.GRCh38.cdna.ncrna_STAR2.5.3a

define startxindexrule
$(1)/SA: $(2)
	mkdir -p $$(@D)
	$(STAR) --runMode genomeGenerate --runThreadN $(nthreads) --genomeDir $(1) \
	--genomeFastaFiles $(2) --limitGenomeGenerateRAM 153509114410
endef
$(eval $(call startxindexrule,$(STARindextxome),$(txome)))
$(foreach F,$(fastqfiles),$(eval $(call startxindexrule,reference/STAR/$(notdir $(F))_stringtie_tx_STAR2.5.3a,stringtie/$(notdir $(F))/$(notdir $(F))_stringtie_tx.fa)))
$(eval $(call startxindexrule,reference/STAR/chess2.0_assembly_fixed_STAR2.5.3a,$(txome_chess)))

## ==================================================================================== ##
##                                 STAR + Salmon                                        ##
## ==================================================================================== ##
## Run STAR to align to the transcriptome
define startxrule
$(4)/$(2)/$(2)_Aligned.out.bam: $(3)/SA $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(STAR) --genomeDir $(3) \
	--readFilesIn $(1)_R1.fastq.gz $(1)_R2.fastq.gz \
	--runThreadN $(nthreads) --outFileNamePrefix $(4)/$(2)/$(2)_ \
	--outFilterMultimapNmax 200 --outFilterMismatchNmax 99999 --outFilterMismatchNoverLmax 0.2 \
	--alignIntronMin 1000 --alignIntronMax 0 --limitOutSAMoneReadBytes 1000000 \
	--outSAMtype BAM Unsorted --readFilesCommand gunzip -c 
endef
$(foreach F,$(fastqfiles),$(eval $(call startxrule,$(F),$(notdir $(F)),$(STARindextxome),STARtxome)))
$(foreach F,$(fastqfiles),$(eval $(call startxrule,$(F),$(notdir $(F)),reference/STAR/$(notdir $(F))_stringtie_tx_STAR2.5.3a,STARtxome_stringtie_tx)))
$(foreach F,$(fastqfilesreal),$(eval $(call startxrule,$(F),$(notdir $(F)),reference/STAR/chess2.0_assembly_fixed_STAR2.5.3a,STARtxome_chess)))

## Run Salmon in alignment-based mode
define salmonstarrule
$(2)/$(1)/quant.sf: $(3)/$(1)/$(1)_Aligned.out.bam $(4)
	mkdir -p $$(@D)
	$(salmon) quant -l A -a $$(word 1,$$^) -t $(4) -p $(nthreads) --seqBias --gcBias -o $$(@D) 
endef
$(foreach F,$(fastqfiles),$(eval $(call salmonstarrule,$(notdir $(F)),salmonstartx,STARtxome,$(txome))))
$(foreach F,$(fastqfiles),$(eval $(call salmonstarrule,$(notdir $(F)),salmonstartx_stringtie_tx,STARtxome_stringtie_tx,stringtie/$(notdir $(F))/$(notdir $(F))_stringtie_tx.fa)))
$(foerach F,$(fastqfilesreal),$(eval $(call salmonstarrule,$(notdir $(F)),salmonstartx_chess,STARtxome_chess,$(txome_chess))))