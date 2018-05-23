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

## Add symbol information to tx2gene
define tx2genesymbolrule
reference/$(1)$(2)_tx2gene_withsymbol.rds: reference/$(1)$(2)_tx2gene.rds $(tx2geneext) Rscripts/add_symbol_to_tx2gene.R
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

