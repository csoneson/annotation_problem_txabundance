bwa := /home/charlotte/software/bwa/bwa
salmon := /home/charlotte/software/Salmon-0.9.1_linux_x86_64/bin/salmon

define bwaindexrule
$(1).sa: $(2)
	mkdir -p $$(@D)
	scp $(2) $(1)
	$(bwa) index $(1)
endef
$(eval $(call bwaindexrule,reference/bwa/Homo_sapiens.GRCh38.cdna.ncrna/Homo_sapiens.GRCh38.cdna.ncrna.fa,$(txome)))
$(foreach F,$(fastqfiles),$(eval $(call bwaindexrule,reference/bwa/$(notdir $(F))_stringtie_tx/$(notdir $(F))_stringtie_tx.fa,stringtie/$(notdir $(F))/$(notdir $(F))_stringtie_tx.fa)))

## ==================================================================================== ##
##                                  BWA + Salmon                                        ##
## ==================================================================================== ##
## Run BWA
define bwarule
$(4)/$(2)/$(2).bam: $(3).sa# $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	$(bwa) mem -t $(nthreads) $(3) $(1)_R1.fastq.gz $(1)_R2.fastq.gz | $(samtools) view -b -F 0x0800 -@ $(nthreads) - > $$@
endef
$(foreach F,$(fastqfiles),$(eval $(call bwarule,$(F),$(notdir $(F)),reference/bwa/Homo_sapiens.GRCh38.cdna.ncrna/Homo_sapiens.GRCh38.cdna.ncrna.fa,bwa/cDNAncRNA)))
$(foreach F,$(fastqfiles),$(eval $(call bwarule,$(F),$(notdir $(F)),reference/bwa/$(notdir $(F))_stringtie_tx/$(notdir $(F))_stringtie_tx.fa,bwa_stringtie_tx)))

## Run Salmon in alignment-based mode
define salmonbwarule
$(2)/$(1)/quant.sf: $(3)/$(1)/$(1).bam $(4)
	mkdir -p $$(@D)
	$(salmon) quant -l A -a $$(word 1,$$^) -t $(4) -p $(nthreads) --seqBias --gcBias -o $$(@D) 
endef
$(foreach F,$(fastqfiles),$(eval $(call salmonbwarule,$(notdir $(F)),salmonbwa/cDNAncRNA,bwa/cDNAncRNA,$(txome))))
$(foreach F,$(fastqfiles),$(eval $(call salmonbwarule,$(notdir $(F)),salmonbwa_stringtie_tx,bwa_stringtie_tx,reference/bwa/$(notdir $(F))_stringtie_tx/$(notdir $(F))_stringtie_tx.fa)))
