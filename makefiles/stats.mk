define covpredsumrule
stats/alpine_coverage_prediction_summary_$(1)$(2).txt: alpine/$(1)$(2)/alpine_predicted_coverage.rds \
Rscripts/alpine_summarize_predicted_coverage.R
	mkdir -p $$(@D)
	$(R) "--args covrds='$$(word 1,$$^)' outtxt='$$@'" Rscripts/alpine_summarize_predicted_coverage.R Rout/alpine_summarize_predicted_coverage_$(1)$(2).Rout
endef
$(foreach F,$(fastqfiles),$(eval $(call covpredsumrule,$(notdir $(F)),)))
$(foreach F,$(fastqfiles),$(eval $(call covpredsumrule,$(notdir $(F)),_stringtie_tx)))
