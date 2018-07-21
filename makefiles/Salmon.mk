salmon := /home/charlotte/software/salmon-0.11.0-linux_x86_64/bin/salmon

salmoncdnancrnaindex := $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.ncrna_sidx_v0.11.0
salmoncdsindex := $(refdir)/cds/Homo_sapiens.GRCh38.cds.all_sidx_v0.11.0
salmonkeepdupindex := reference/salmon/Homo_sapiens.TRCh38.cdna.ncrna.keepdup_sidx_v0.11.0

## Build Salmon index
define salmonindexrule
$(1)/hash.bin: $(3)
	$(2) index -t $$< -k 25 -i $$(@D) -p $(nthreads) --type quasi $(4)
endef
$(eval $(call salmonindexrule,$(salmoncdnancrnaindex),$(salmon),$(txome),))
$(eval $(call salmonindexrule,$(salmoncdsindex),$(salmon),$(cds),))
$(eval $(call salmonindexrule,$(salmonkeepdupindex),$(salmon),$(txome),--keepDuplicates))
$(foreach F,$(fastqfiles),$(eval $(call salmonindexrule,reference/salmon/$(notdir $(F))_stringtie_tx_sidx_v0.11.0,$(salmon),stringtie/$(notdir $(F))/$(notdir $(F))_stringtie_tx.fa,)))
$(eval $(call salmonindexrule,reference/salmon/chess2.0_assembly_fixed_sidx_v0.11.0,$(salmon),$(txome_chess),))
$(eval $(call salmonindexrule,reference/salmon/chess2.0_assembly_fixed_keepdup_sidx_v0.11.0,$(salmon),$(txome_chess),--keepDuplicates))

## Count number of removed transcripts
stats/nbr_duplicate_transcripts_Salmon_Homo_sapiens.GRCh38.cdna.ncrna_sidx_v0.11.0.txt: $(salmoncdnancrnaindex)/hash.bin Rscripts/count_salmon_duplicate_tx.R
	$(R) "--args salmonindexdir='$(salmoncdnancrnaindex)' tx2gene='$(tx2geneext)' outtxt='$@'" Rscripts/count_salmon_duplicate_tx.R Rout/count_salmon_duplicate_tx_Homo_sapiens.GRCh38.cdna.ncrna_sidx_v0.11.0.Rout

## Run Salmon
define salmonrule
$(5)/$(2)/quant.sf: $(3)/hash.bin $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(4) quant -i $$(word 1,$$(^D)) -l A -p $(nthreads) -1 $(1)_R1.fastq.gz -2 $(1)_R2.fastq.gz -o $$(@D) --seqBias --gcBias --posBias $(6)
endef
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(notdir $(F)),$(salmoncdnancrnaindex),$(salmon),salmon/cDNAncRNA,--numBootstraps 100)))
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(notdir $(F)),$(salmoncdsindex),$(salmon),salmon/cds,)))
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(notdir $(F)),$(salmonkeepdupindex),$(salmon),salmon/cDNAncRNAkeepdup,)))
$(foreach F,$(fastqfiles),$(eval $(call salmonrule,$(F),$(notdir $(F)),reference/salmon/$(notdir $(F))_stringtie_tx_sidx_v0.11.0,$(salmon),salmon_stringtie_tx,)))
$(foreach F,$(fastqfilesreal),$(eval $(call salmonrule,$(F),$(notdir $(F)),reference/salmon/chess2.0_assembly_fixed_sidx_v0.11.0,$(salmon),salmon_chess,--numBootstraps 100)))
$(foreach F,$(fastqfilesreal),$(eval $(call salmonrule,$(F),$(notdir $(F)),reference/salmon/chess2.0_assembly_fixed_keepdup_sidx_v0.11.0,$(salmon),salmon_chesskeepdup,)))