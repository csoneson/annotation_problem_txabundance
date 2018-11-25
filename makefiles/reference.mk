## Reference files
refdir := /home/Shared/data/annotation/Human/Ensembl_GRCh38.90
cdna := $(refdir)/cDNA/Homo_sapiens.GRCh38.cdna.all.fa
cds := $(refdir)/cds/Homo_sapiens.GRCh38.cds.all.fa
ncrna := $(refdir)/ncRNA/Homo_sapiens.GRCh38.ncrna.fa
genome := $(refdir)/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf := $(refdir)/gtf/Homo_sapiens.GRCh38.90.gtf
txome := reference/Homo_sapiens.GRCh38.cdna.ncrna.fa

gvizgenemodels := reference/Gviz/Homo_sapiens.GRCh38.90_gviz_genemodels.rds

tx2gene := reference/Homo_sapiens.GRCh38.90_tx2gene.rds
tx2geneext := reference/Homo_sapiens.GRCh38.90_tx2gene_ext.rds
flatgtfexons := reference/Homo_sapiens.GRCh38.90.reduced.exons.gtf
flatgtfintrons := reference/Homo_sapiens.GRCh38.90.introns.gtf

gffread := /home/charlotte/software/gffread-0.9.12.Linux_x86_64/gffread

$(txome): $(cdna) $(ncrna)
	cat $(cdna) $(ncrna) > $(txome)

## Generate tx2gene
$(tx2gene): $(cdna) $(cds) $(ncrna) Rscripts/generate_tx2gene.R
	$(R) "--args cdna='$(word 1,$^)' cds='$(word 2,$^)' ncrna='$(word 3,$^)' outrds='$@'" Rscripts/generate_tx2gene.R Rout/generate_tx2gene.Rout 

$(tx2geneext): $(tx2gene) $(gtf) Rscripts/extend_tx2gene.R
	$(R) "--args tx2gene='$(tx2gene)' gtf='$(gtf)' outrds='$@'" Rscripts/extend_tx2gene.R Rout/extend_tx2gene.Rout

## Flatten the gtf file and generate a separate gtf file with introns (gene range\union of exons)
$(flatgtfintrons): $(gtf) Rscripts/generate_exon_and_intron_gtfs.R
	$(R) "--args ingtf='$(gtf)' outexongtf='$(flatgtfexons)' outintrongtf='$(flatgtfintrons)'" Rscripts/generate_exon_and_intron_gtfs.R Rout/generate_exon_and_intron_gtfs.Rout

$(flatgtfexons): $(flatgtfintrons)
	touch $(flatgtfexons)

## ==================================================================================== ##
##                     Reference files for extended annotation                          ##
## ==================================================================================== ##
## Generate tx2gene
define tx2generule
reference/$(1)$(2)_tx2gene.rds: $(3) Rscripts/generate_tx2gene_from_gtf.R
	$(R) "--args gtf='$$<' outrds='$$@'" Rscripts/generate_tx2gene_from_gtf.R Rout/generate_tx2gene_from_gtf_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call tx2generule,$(notdir $(F)),_stringtie_tx,stringtie/$(notdir $(F))/$(notdir $(F))_filtered.gtf)))

define gene2symbolrule
reference/$(1)$(2)_gene2symbol.rds: $(3) Rscripts/generate_gene2symbol_from_gtf.R
	$(R) "--args gtf='$$<' outrds='$$@'" Rscripts/generate_gene2symbol_from_gtf.R Rout/generate_gene2symbol_from_gtf_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call gene2symbolrule,$(notdir $(F)),_stringtie_tx,stringtie/$(notdir $(F))/$(notdir $(F))_filtered.gtf)))

define tx2symbolrule
reference/$(1)$(2)_tx2symbol.rds: $(3) Rscripts/generate_tx2symbol_from_gtf.R
	$(R) "--args gtf='$$<' outrds='$$@'" Rscripts/generate_tx2symbol_from_gtf.R Rout/generate_tx2symbol_from_gtf_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call tx2symbolrule,$(notdir $(F)),_stringtie_tx,stringtie/$(notdir $(F))/$(notdir $(F))_filtered.gtf)))

## Add symbol information to tx2gene
define tx2genesymbolrule
reference/$(1)$(2)_tx2gene_withsymbol.rds: reference/$(1)$(2)_tx2gene.rds reference/$(1)$(2)_gene2symbol.rds Rscripts/add_symbol_to_tx2gene.R
	$(R) "--args tx2gene='$$(word 1,$$^)' info='$$(word 2,$$^)' outrds='$$@'" Rscripts/add_symbol_to_tx2gene.R Rout/add_symbol_to_tx2gene_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call tx2genesymbolrule,$(notdir $(F)),_stringtie_tx)))

## ==================================================================================== ##
##                            characterize genes                                        ##
## ==================================================================================== ##
output/gene_characteristics.rds: $(gtf) $(txome) Rscripts/calculate_gene_characteristics.R
	mkdir -p $(@D)
	$(R) "--args gtf='$(gtf)' txome='$(txome)' outrds='$@'" Rscripts/calculate_gene_characteristics.R Rout/calculate_gene_characteristics.Rout

## ==================================================================================== ##
##                              Gviz gene models                                        ##
## ==================================================================================== ##
## Gene models for Gviz
$(gvizgenemodels): $(gtf) Rscripts/generate_genemodels.R Rscripts/helper_plot_tracks.R
	mkdir -p $(@D)
	$(R) "--args gtf='$(gtf)' outrds='$@'" Rscripts/generate_genemodels.R Rout/generate_genemodels.Rout

## ==================================================================================== ##
##                                CHESS annotation                                      ##
## ==================================================================================== ##
gtf_chess := reference_chess/chess2.0_assembly_fixed.gtf
txome_chess := reference_chess/chess2.0_assembly_fixed.fa
tx2gene_chess := reference_chess/chess2.0_assembly_fixed_tx2gene.rds
tx2gene_chess_withsymbol := reference_chess/chess2.0_assembly_fixed_tx2gene_withsymbol.rds
info_chess := reference_chess/gene_id_to_symbol.rds
gvizgenemodels_chess := reference_chess/chess2.0_assembly_fixed_gviz_genemodels.rds

$(gtf_chess): reference_chess/chess2.0_assembly.gff reference_chess/chess2.0.genes Rscripts/fix_chess_gtf.R
	$(R) "--args ingtf='$(word 1,$^)' ingenes='$(word 2,$^)' outgtf='$@' outinfo='$(info_chess)'" Rscripts/fix_chess_gtf.R Rout/fix_chess_gtf.Rout

$(info_chess): $(gtf_chess)
	touch $@

$(txome_chess): $(gtf_chess) $(genome)
	$(gffread) -w $@ -g $(genome) $<

$(tx2gene_chess): $(gtf_chess) Rscripts/generate_tx2gene_from_gtf.R
	$(R) "--args gtf='$<' outrds='$@'" Rscripts/generate_tx2gene_from_gtf.R Rout/generate_tx2gene_from_gtf_chess.Rout

$(tx2gene_chess_withsymbol): $(tx2gene_chess) $(info_chess) Rscripts/add_symbol_to_tx2gene.R
	$(R) "--args tx2gene='$(tx2gene_chess)' info='$(info_chess)' outrds='$@'" Rscripts/add_symbol_to_tx2gene.R Rout/add_symbol_to_tx2gene_chess.Rout

$(gvizgenemodels_chess): $(gtf_chess) Rscripts/generate_genemodels.R Rscripts/helper_plot_tracks.R
	mkdir -p $(@D)
	$(R) "--args gtf='$(gtf_chess)' outrds='$@'" Rscripts/generate_genemodels.R Rout/generate_genemodels_chess.Rout

output/gene_characteristics_chess.rds: $(gtf_chess) Rscripts/calculate_gene_characteristics_chess.R
	mkdir -p $(@D)
	$(R) "--args gtf='$(gtf_chess)' outrds='$@'" Rscripts/calculate_gene_characteristics_chess.R Rout/calculate_gene_characteristics_chess.Rout


## ==================================================================================== ##
##                    Replace shorter 3'UTRs with longer ones                           ##
## ==================================================================================== ##
gtf_longutr_added := reference_longUTR_added/Homo_sapiens.GRCh38.90_longUTR_added.gtf
txome_longutr_added := reference_longUTR_added/Homo_sapiens.GRCh38.90_longUTR_added.fasta
tx2gene_longutr_added := reference_longUTR_added/Homo_sapiens.GRCh38.90_longUTR_added_tx2gene.rds
gvizgenemodels_longutr_added := reference_longUTR_added/Homo_sapiens.GRCh38.90_longUTR_added_gviz_genemodels.rds

$(txome_longutr_added): $(gtf) $(txome) Rscripts/extend_annotation_longest_3putr.R
	mkdir -p $(@D)
	$(R) "--args gtffile='$(gtf)' txfastafile='$(txome)' outbase='reference_longUTR_added/Homo_sapiens.GRCh38.90_longUTR_added'" Rscripts/extend_annotation_longest_3putr.R Rout/extend_annotation_longest_3putr.Rout

$(gtf_longutr_added): $(txome_longutr_added)
	touch $@

$(tx2gene_longutr_added): $(txome_longutr_added)
	touch $@

$(gvizgenemodels_longutr_added): $(gtf_longutr_added) Rscripts/generate_genemodels.R Rscripts/helper_plot_tracks.R
	mkdir -p $(@D)
	$(R) "--args gtf='$(gtf_longutr_added)' outrds='$@'" Rscripts/generate_genemodels.R Rout/generate_genemodels_longUTR_added.Rout

