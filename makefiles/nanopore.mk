minimap2 := /home/charlotte/software/minimap2_2.10-r773-dirty/minimap2
bam_count_reads := /home/charlotte/miniconda3/bin/bam_count_reads.py

nanopore_dir := /home/Shared/data/seq/hussain_bath_nanopore_rnaseq/FGCZ/FASTQ
nanopore_sample := 20171207_1645_p2557_4017_2_ALLREADS.pass

## Align to transcriptome
minimap2txome/$(nanopore_sample)/$(nanopore_sample)_minimap2_txome.bam: $(nanopore_dir)/$(nanopore_sample).FASTQ.gz $(txome)
	mkdir -p $(@D)
	$(minimap2) -t $(nthreads) -ax map-ont $(txome) $< | $(samtools) view -bS - > $@

## Apply wub to get transcript abundance estimates
minimap2wub/$(nanopore_sample)/bam_count_reads.tsv: minimap2txome/$(nanopore_sample)/$(nanopore_sample)_minimap2_txome.bam $(txome)
	mkdir -p $(@D)
	$(bam_count_reads) -a 5 -z $(txome) -t $@ $<

## Apply Salmon to get transcript abundance estimates
minimap2salmon/$(nanopore_sample)/quant.sf: minimap2txome/$(nanopore_sample)/$(nanopore_sample)_minimap2_txome.bam $(txome)
	mkdir -p $(@D)
	$(salmon) quant -t $(txome) -l A -a $< -o $(@D) -p $(nthreads) --fldMax 230000 --fldMean 600

## Align to genome
minimap2genome/$(nanopore_sample)/$(nanopore_sample)_minimap2_genome_s.bam: $(nanopore_dir)/$(nanopore_sample).FASTQ.gz $(genome)# $(minimap2)
	mkdir -p $(@D)
	$(minimap2) -t $(nthreads) -ax splice $(genome) $< | $(samtools) view -bS - | \
	$(samtools) sort -o $@ -T tmp$(1) -@ $(nthreads) - 

## Index genome bam
minimap2genome/$(nanopore_sample)/$(nanopore_sample)_minimap2_genome_s.bam.bai: minimap2genome/$(nanopore_sample)/$(nanopore_sample)_minimap2_genome_s.bam
	$(samtools) index $<

## Convert BAM files to bigWig
minimap2genomebigwig/$(nanopore_sample)_minimap2_genome_s.bw: $(STARindexnogtf)/chrNameLength.txt \
minimap2genome/$(nanopore_sample)/$(nanopore_sample)_minimap2_genome_s.bam
	mkdir -p $(@D)	
	$(bedtools) genomecov -split -ibam $(word 2,$^) \
	-bg > minimap2genomebigwig/$(nanopore_sample)_minimap2_genome_s.bedGraph
	
	bedGraphToBigWig minimap2genomebigwig/$(nanopore_sample)_minimap2_genome_s.bedGraph \
	$(word 1,$^) $@
	
	rm -f minimap2genomebigwig/$(nanopore_sample)_minimap2_genome_s.bedGraph

