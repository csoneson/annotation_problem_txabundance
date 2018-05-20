## Define new panel function for ggpairs, to show both Pearson and Spearman correlation
combinecor <- function(data, mapping, color = I("black"), sizeRange = c(1, 5), ...) {
  
  # get the x and y data to use the other code
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  
  ct1 <- cor(x, y, method = "pearson", use = "pairwise.complete.obs")
  ct2 <- cor(x, y, method = "spearman", use = "pairwise.complete.obs")
  
  # since we can't print it to get the strsize, just use the max size range
  cex <- max(sizeRange)
  
  
  # plot the cor value
  ggally_text(
    label = paste0("Pearson: \n", signif(ct1, 3), "\nSpearman: \n", signif(ct2, 3)), 
    mapping = mapping,
    xP = 0.5, yP = 0.5, 
    color = color,
    ...
  ) 
}
