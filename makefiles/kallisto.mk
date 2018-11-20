kallisto := /home/charlotte/software/kallisto_linux-v0.44.0/kallisto

kallistocdnancrnaindex := $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.ncrna.all_kidx_v0.44.0

## Build kallisto index
define kallistoindexrule
$(1): $(2)
	mkdir -p $$(@D)
	$(kallisto) index -i $$@ -k 25 $$<
endef
$(eval $(call kallistoindexrule,$(kallistocdnancrnaindex),$(txome)))
$(foreach F,$(fastqfiles),$(eval $(call kallistoindexrule,reference/kallisto/$(notdir $(F))_stringtie_tx_kidx_v0.44.0,stringtie/$(notdir $(F))/$(notdir $(F))_stringtie_tx.fa)))
$(eval $(call kallistoindexrule,reference/kallisto/chess2.0_assembly_fixed_kidx_v0.44.0,$(txome_chess)))
$(eval $(call kallistoindexrule,reference/kallisto/Homo_sapiens.GRCh38.cdna.ncrna_longUTR_added_kidx_v0.44.0,$(txome_longutr_added)))

## ==================================================================================== ##
##                                   kallisto                                           ##
## ==================================================================================== ##
## Run kallisto
define kallistorule
$(4)/$(2)/abundance.tsv: $(3) $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(kallisto) quant -i $(3) -o $$(@D) --bias -t $(nthreads) $(5) $(1)_R1.fastq.gz $(1)_R2.fastq.gz
endef
$(foreach F,$(fastqfiles),$(eval $(call kallistorule,$(F),$(notdir $(F)),$(kallistocdnancrnaindex),kallisto/cDNAncRNA,--rf-stranded)))
$(foreach F,$(fastqfiles),$(eval $(call kallistorule,$(F),$(notdir $(F)),reference/kallisto/$(notdir $(F))_stringtie_tx_kidx_v0.44.0,kallisto_stringtie_tx,--rf-stranded)))
$(foreach F,$(fastqfilesreal),$(eval $(call kallistorule,$(F),$(notdir $(F)),reference/kallisto/chess2.0_assembly_fixed_kidx_v0.44.0,kallisto_chess,--rf-stranded)))
$(foreach F,$(fastqfilesreal),$(eval $(call kallistorule,$(F),$(notdir $(F)),reference/kallisto/Homo_sapiens.GRCh38.cdna.ncrna_longUTR_added_kidx_v0.44.0,kallisto_longUTR_added,--rf-stranded)))
