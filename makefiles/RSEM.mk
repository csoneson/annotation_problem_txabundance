RSEM := /home/charlotte/software/RSEM-1.3.0

define gene2txrule
reference/$(1)_rsemgene2tx.txt: $(2) Rscripts/generate_rsemgene2tx.R
	$(R) "--args tx2gene='$(2)' rsemgene2tx='$$@'" Rscripts/generate_rsemgene2tx.R Rout/generate_rsemgene2tx_$(1).Rout
endef
$(eval $(call gene2txrule,Homo_sapiens.GRCh38.90,$(tx2gene)))
$(foreach F,$(fastqfiles),$(eval $(call gene2txrule,$(notdir $(F))_stringtie_tx,reference/$(notdir $(F))_stringtie_tx_tx2gene.rds)))

define rsemindexrule
reference/RSEM/$(1)/$(1).n2g.idx.fa: $(2) $(3)
	mkdir -p $$(@D)
	$(RSEM)/rsem-prepare-reference --transcript-to-gene-map $(3) --bowtie --bowtie-path /usr/bin \
	$(2) reference/RSEM/$(1)/$(1)
endef
$(eval $(call rsemindexrule,Homo_sapiens.GRCh38.rsem.cdna.ncrna,$(txome),reference/Homo_sapiens.GRCh38.90_rsemgene2tx.txt))
$(foreach F,$(fastqfiles),$(eval $(call rsemindexrule,$(notdir $(F))_stringtie_tx,stringtie/$(notdir $(F))/$(notdir $(F))_stringtie_tx.fa,reference/$(notdir $(F))_stringtie_tx_rsemgene2tx.txt)))

## ==================================================================================== ##
##                                     RSEM                                             ##
## ==================================================================================== ##
## Run RSEM
define rsemrule
$(4)/$(2)/$(2).isoforms.results: $(3).n2g.idx.fa# $(1)_R1.fastq.gz $(1)_R2.fastq.gz
	mkdir -p $$(@D)
	bash -c '$(RSEM)/rsem-calculate-expression -p $(nthreads) --bowtie-path /usr/bin --paired-end <(gunzip -c $(1)_R1.fastq.gz) <(gunzip -c $(1)_R2.fastq.gz) $(3) $(4)/$(2)/$(2)'
	rm -f $(4)/$(2)/$(2).transcript.bam
endef
$(foreach F,$(fastqfiles),$(eval $(call rsemrule,$(F),$(notdir $(F)),reference/RSEM/Homo_sapiens.GRCh38.rsem.cdna.ncrna/Homo_sapiens.GRCh38.rsem.cdna.ncrna,RSEM/cDNAncRNA)))
$(foreach F,$(fastqfiles),$(eval $(call rsemrule,$(F),$(notdir $(F)),reference/RSEM/$(notdir $(F))_stringtie_tx/$(notdir $(F))_stringtie_tx,RSEM_stringtie_tx)))
