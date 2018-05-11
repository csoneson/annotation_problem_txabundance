simulation/misannotated_utr/sim_misannotated_utr_1_R2.fastq.gz: $(gtf) $(txome) Rscripts/simulate_reads_with_misannotated_utr.R
	mkdir -p $(@D)
	$(R) "--args gtf='$(gtf)' txfasta='$(txome)' outfasta='$(@D)/modified_transcripts.fa' readdir='$(@D)' readbasename='sim_misannotated_utr_1'" Rscripts/simulate_reads_with_misannotated_utr.R Rout/simulate_reads_with_misannotated_utr.Rout

simulation/misannotated_utr/sim_misannotated_utr_1_R1.fastq.gz: simulation/misannotated_utr/sim_misannotated_utr_1_R2.fastq.gz
	touch $@