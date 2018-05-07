hera := /home/charlotte/software/hera/build/hera
hera_build := /home/charlotte/software/hera/build/hera_build

define heraindexrule
reference/hera/$(1)/index: $(genome) $(2)
	mkdir -p $$(@D)
	$(hera_build) --fasta $(genome) --gtf $(2) --outdir $$(@D)/
endef
$(eval $(call heraindexrule,Homo_sapiens.GRCh38,$(gtf)))
$(foreach F,$(fastqfiles),$(eval $(call heraindexrule,$(notdir $(F))_stringtie_tx,stringtie/$(notdir $(F))/$(notdir $(F))_filtered_withgene.gtf)))

## ==================================================================================== ##
##                                      HERA                                            ##
## ==================================================================================== ##
## Run hera
define herarule
hera$(3)/$(2)/abundance.tsv: $(4)/index# $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(hera) quant -i $(4) -o $$(@D) -w 1 -t $(nthreads) $(1)_R1.fastq.gz $(1)_R2.fastq.gz
endef
$(foreach F,$(fastqfiles),$(eval $(call herarule,$(F),$(notdir $(F)),,reference/hera/Homo_sapiens.GRCh38)))
$(foreach F,$(fastqfiles),$(eval $(call herarule,$(F),$(notdir $(F)),_stringtie_tx,reference/hera/$(notdir $(F))_stringtie_tx)))
