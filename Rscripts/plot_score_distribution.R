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

png(gsub("rds$", "png", outrds), width = 14, 
    height = 14, unit = "in", res = 300)
  
## Plot score distribution
print(
  cowplot::plot_grid(
    ggplot(scores, aes(x = 1, y = score, fill = method)) + 
      geom_violin() + theme_bw() + facet_wrap(~ method, nrow = 1) + 
      scale_color_manual(values = method_colors) + 
      theme(legend.position = "none") + xlab(""), 
    
    ggplot(scores, aes(x = uniqjuncfraction, y = score, color = method)) +
      geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) +
      geom_smooth(color = "black", method = "loess") + xlab("Fraction uniquely mapping junction reads") +
      theme(legend.position = "none"),
    
    ggplot(scores %>% dplyr::mutate(intron_exon_ratio = replace(intron_exon_ratio, 
                                                                intron_exon_ratio > 10, 10)), 
           aes(x = intron_exon_ratio, y = score, color = method)) + 
      geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) + 
      geom_smooth(color = "black", method = "loess") + xlab("Intron/exon count ratio") + 
      theme(legend.position = "none"),
      
    ggplot(scores, aes(x = count, y = score, color = method)) + 
      geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) + 
      geom_smooth(color = "black", method = "loess") + xlab("Gene count") + scale_x_sqrt() + 
      theme(legend.position = "none"),
      
    ggplot(scores, aes(x = length_diff_3putrs_samestart, y = score, color = method)) + 
      geom_point(alpha = 0.5, size = 0.5) + theme_bw() + facet_wrap(~ method, nrow = 1) + 
      geom_smooth(color = "black", method = "loess") + xlab("Length difference of 3'UTRs with same start") + 
      theme(legend.position = "none"),
    
    ncol = 1      
  )
)

dev.off()

saveRDS(NULL, file = outrds)
sessionInfo()
date()
