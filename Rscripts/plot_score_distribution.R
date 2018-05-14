args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

quantmethods <- strsplit(quantmethods, ",")[[1]]

print(scorerds)
print(quantmethods)
print(outrds)  ## output file

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

scores <- readRDS(scorerds)$genes

## Define colors
method_colors <- c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D",
                   "#4EB265", "#CAEDAB", "#777777", "#E8601C",
                   "#1965B0", "#882E72", "#F6C141", "#F7EE55",
                   "#90C987")[seq_len(length(unique(scores$method)))]
names(method_colors) <- unique(scores$method)

scores <- scores %>% dplyr::filter(method %in% quantmethods)

## Characterize genes with NA score
scoresNA <- scores %>% #dplyr::filter(is.na(score)) %>% 
  dplyr::group_by(method) %>% 
  dplyr::summarize(
    valid_score = sum(!is.na(score)),
    not_expressed = sum(is.na(score) & count == 0),
    no_junctions = sum(is.na(score) & count > 0 & nbr_junctions_in_gene == 0, na.rm = TRUE),
    junctions_but_no_unique_reads = sum(is.na(score) & count > 0 & 
                                          nbr_junctions_in_gene > 0 & uniqjuncreads == 0, na.rm = TRUE),
    junctions_but_too_few_unique_reads = sum(is.na(score) & count > 0 & nbr_junctions_in_gene > 0 & 
                                               uniqjuncreads > 0 & uniqjuncfraction <= 0.75, na.rm = TRUE),
    other = length(count) - valid_score - not_expressed - no_junctions - junctions_but_no_unique_reads - 
      junctions_but_too_few_unique_reads) %>% 
  tidyr::gather(reason, number, -method) %>% 
  dplyr::mutate(reason = replace(reason, reason == "valid_score", "Valid score"),
                reason = replace(reason, reason == "not_expressed", "Estimated abundance = 0"),
                reason = replace(reason, reason == "no_junctions", "No junctions"),
                reason = replace(reason, reason == "junctions_but_no_unique_reads", 
                                 "Has junctions, but no uniquely mapping junction reads"),
                reason = replace(reason, reason == "junctions_but_too_few_unique_reads",
                                 "Has junctions, but too large fraction multimapping junction reads"),
                reason = replace(reason, reason == "other", "Other")) %>%
  dplyr::mutate(reason = factor(reason, levels = c("Has junctions, but too large fraction multimapping junction reads",
                                                   "Has junctions, but no uniquely mapping junction reads",
                                                   "No junctions",
                                                   "Other", "Estimated abundance = 0", "Valid score")))
png(gsub("\\.rds$", "_NAscores.png", outrds), width = 7, height = 6, 
    unit = "in", res = 300)
print(ggplot(scoresNA, aes(x = method, y = number, fill = reason)) + 
        geom_bar(stat = "identity") + theme_bw() + 
        theme(legend.position = "bottom") + 
        guides(fill = guide_legend(nrow = 2, title = "")) + 
        scale_fill_manual(values = c("blue", "orange", "green", "cyan", "black", "grey")) + 
        xlab("") + ylab("Number of genes"))
dev.off()

for (mincount in c(0, 1000)) {
  png(gsub("\\.rds$", paste0("_min", mincount, "reads.png"), outrds), width = 14, 
      height = 14, unit = "in", res = 300)
  
  ## Plot score distribution
  print(
    cowplot::plot_grid(
      ggplot(scores %>% dplyr::filter(!is.na(score) & count >= mincount), 
             aes(x = 1, y = score, fill = method)) + 
        geom_violin() + theme_bw() + facet_wrap(~ method, nrow = 1) + 
        scale_fill_manual(values = method_colors) + 
        theme(legend.position = "none") + xlab(""), 
      
      ggplot(scores %>% dplyr::filter(!is.na(score) & count >= mincount),
             aes(x = uniqjuncfraction, y = score, color = method)) +
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) +
        geom_smooth(color = "black", method = "loess") + xlab("Fraction uniquely mapping junction reads") +
        theme(legend.position = "none") + 
        scale_color_manual(values = method_colors),
      
      ggplot(scores %>% dplyr::filter(!is.na(score) & count >= mincount),
             aes(x = uniqjuncreads, y = score, color = method)) +
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) +
        geom_smooth(color = "black", method = "loess") + xlab("Number of uniquely mapping junction reads") +
        theme(legend.position = "none") + 
        scale_color_manual(values = method_colors),
      
      ggplot(scores %>% dplyr::filter(!is.na(score) & count >= mincount) %>% 
               dplyr::mutate(intron_exon_ratio = replace(intron_exon_ratio, 
                                                         intron_exon_ratio > 10, 10)), 
             aes(x = intron_exon_ratio, y = score, color = method)) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) + 
        geom_smooth(color = "black", method = "loess") + xlab("Intron/exon count ratio") + 
        theme(legend.position = "none") + 
        scale_color_manual(values = method_colors),
      
      ggplot(scores %>% dplyr::filter(!is.na(score) & count >= mincount),
             aes(x = count, y = score, color = method)) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1, scales = "free_x") + 
        geom_smooth(color = "black", method = "loess") + xlab("Gene count") + scale_x_sqrt() + 
        theme(legend.position = "none") + 
        scale_color_manual(values = method_colors),
    
      ggplot(scores %>% dplyr::filter(!is.na(score) & count >= mincount) %>%
               dplyr::mutate(nbr_transcripts = replace(nbr_transcripts,
                                                       nbr_transcripts > 60, 60)), 
             aes(x = nbr_transcripts, y = score, color = method)) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) + 
        geom_smooth(color = "black", method = "loess") + xlab("Number of transcripts") +  
        theme(legend.position = "none") + 
        scale_color_manual(values = method_colors),
      
      ggplot(scores %>% dplyr::filter(!is.na(score) & count >= mincount) %>%
               dplyr::mutate(nbr_junctions_in_gene = replace(nbr_junctions_in_gene,
                                                             nbr_junctions_in_gene > 700, 700)),
             aes(x = nbr_junctions_in_gene, y = score, color = method)) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) + 
        geom_smooth(color = "black", method = "loess") + xlab("Number of junctions") +  
        theme(legend.position = "none") + 
        scale_color_manual(values = method_colors),
      
      ggplot(scores %>% dplyr::filter(!is.na(score) & count >= mincount),
             aes(x = length_diff_3putrs_samestart, y = score, color = method)) + 
        geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) + 
        geom_smooth(color = "black", method = "loess") + xlab("Length difference of 3'UTRs with same start") + 
        theme(legend.position = "none") + 
        scale_color_manual(values = method_colors),
      
      ncol = 1      
    )
  )
  
  dev.off()
}

saveRDS(NULL, file = outrds)
sessionInfo()
date()
