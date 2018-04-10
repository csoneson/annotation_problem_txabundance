salmon := /home/charlotte/software/Salmon-0.9.1_linux_x86_64/bin/salmon

salmoncdnancrnaindex := $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.ncrna_sidx_v0.9.1
salmoncdsindex := $(refdir)/cds/Homo_sapiens.GRCh38.cds.all_sidx_v0.9.1

## Build Salmon index
define salmonindexrule
$(1)/hash.bin: $(2)
	$(salmon) index -t $$< -k 25 -i $$(@D) -p $(nthreads) --type quasi
endef
$(eval $(call salmonindexrule,$(salmoncdnancrnaindex),$(txome)))
$(eval $(call salmonindexrule,$(salmoncdsindex),$(cds)))
$(foreach F,$(fastqfiles),$(eval $(call salmonindexrule,reference/salmon/$(notdir $(F))_stringtie_tx_sidx_v0.9.1,stringtie/$(notdir $(F))/$(notdir $(F))_stringtie_tx.fa)))

## Run Salmon
define salmonrule
$(4)/$(2)/quant.sf: $(3)/hash.bin# $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(salmon) quant -i $$(word 1,$$(^D)) -l A -p $(nthreads) -1 $(1)_R1.fastq.gz -2 $(1)_R2.fastq.gz -o $$(@D) --seqBias --gcBias --posBias
endef
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(notdir $(F)),$(salmoncdnancrnaindex),salmon/cDNAncRNA)))
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(notdir $(F)),$(salmoncdsindex),salmon/cds)))
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(notdir $(F)),reference/salmon/$(notdir $(F))_stringtie_tx_sidx_v0.9.1,salmon_stringtie_tx)))

