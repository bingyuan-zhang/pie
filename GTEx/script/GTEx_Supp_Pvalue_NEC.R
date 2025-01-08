enrich_res <- readRDS("./GTEx/results/enriched_terms.rds")

library(ggplot2)

topN_values <- seq(10, 400, by=10)
overlaps <- c()
for (i in seq_along(topN_values)) {
  n <- topN_values[i]
  set1 <- order(enrich_res$pval)[1:n]
  set2 <- order(enrich_res$nes, decreasing=TRUE)[1:n]
  overlaps[i] <- length(intersect(set1, set2))
}


df <- data.frame(topN_values = topN_values, overlaps = overlaps)

ggplot(df, aes(x = topN_values, y = overlaps)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Top N",
    y = "Number of overlapping terms",
    title = "Overlap in top terms by p-value and NES"
  ) +
  annotate("segment",
           x = 10, xend = 400,
           y = 10, yend = 400,
           colour = "red") +
  annotate("text",
           x = 200,
           y = 250,
           label = "y = x",
           colour = "red",
           size = 5) +
  theme_classic()
