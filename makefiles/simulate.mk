## Simulate reads
simulation/misannotated_utr/sim_misannotated_utr_1_modified_genes.rds: $(gtf) $(txome) Rscripts/simulate_reads_with_misannotated_utr.R
	mkdir -p $(@D)
	$(R) "--args gtf='$(gtf)' txfasta='$(txome)' outfasta='$(@D)/modified_transcripts.fa' readdir='$(@D)' readbasename='sim_misannotated_utr_1' readlen=125" Rscripts/simulate_reads_with_misannotated_utr.R Rout/simulate_reads_with_misannotated_utr.Rout

simulation/misannotated_utr/sample_01_1.fasta.gz: simulation/misannotated_utr/sim_misannotated_utr_1_modified_genes.rds
	touch $@

simulation/misannotated_utr/sample_01_2.fasta.gz: simulation/misannotated_utr/sim_misannotated_utr_1_modified_genes.rds
	touch $@

## Convert fasta files to fastq
simulation/misannotated_utr/sim_misannotated_utr_1_R2.fastq.gz: simulation/misannotated_utr/sample_01_1.fasta.gz \
simulation/misannotated_utr/sample_01_2.fasta.gz
	## Swap first and second reads, to generate strandedness consistent with the other data sets
	zcat simulation/misannotated_utr/sample_01_1.fasta.gz | perl fasta_to_fastq.pl - > simulation/misannotated_utr/sim_misannotated_utr_1_R2_tmp.fastq
	zcat simulation/misannotated_utr/sample_01_2.fasta.gz | perl fasta_to_fastq.pl - > simulation/misannotated_utr/sim_misannotated_utr_1_R1_tmp.fastq

	bbmap/shuffle.sh in=simulation/misannotated_utr/sim_misannotated_utr_1_R1_tmp.fastq \
	in2=simulation/misannotated_utr/sim_misannotated_utr_1_R2_tmp.fastq \
	out=simulation/misannotated_utr/sim_misannotated_utr_1_R1.fastq \
	out2=simulation/misannotated_utr/sim_misannotated_utr_1_R2.fastq
	rm -f simulation/misannotated_utr/sim_misannotated_utr_1_R1_tmp.fastq
	rm -f simulation/misannotated_utr/sim_misannotated_utr_1_R2_tmp.fastq

	gzip simulation/misannotated_utr/sim_misannotated_utr_1_R1.fastq
	gzip simulation/misannotated_utr/sim_misannotated_utr_1_R2.fastq

simulation/misannotated_utr/sim_misannotated_utr_1_R1.fastq.gz: simulation/misannotated_utr/sim_misannotated_utr_1_R2.fastq.gz
	touch $@
