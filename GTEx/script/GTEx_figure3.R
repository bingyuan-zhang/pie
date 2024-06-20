library(ggplot2)
library(dplyr)

emb    <- readRDS("./GTEx/results/GTEx_emb.rds")
terms  <- readRDS("./GTEx/results/enriched_terms.rds")
term_by_intv <- readRDS("./GTEx/results//term_intv.rds")

## Enriched GO terms
terms_pd <- terms$pd
terms_pval <- terms$pval
terms_egene <- terms$egene
terms_name <- terms$name %>%
  strsplit(split = "\\(") %>%
  map_vec(~ head(.x, 1)) %>%
  str_remove("\\)") %>%
  str_remove(" $")

## show the following terms
show_name <-
  c("Striated Muscle Contraction",
    "Signal Release From Synapse",
    "Oxygen Transport",
    "Negative Regulation Of Blood Coagulation",
    "Regulation Of Calcium Ion Transmembrane Transport",
    "Epidermis Development",
    "Extracellular Matrix Organization",
    "Negative Regulation Of Smooth Muscle Cell Proliferation"
  )

show_label <-
  c("Muscle",
    "Brain",
    "Whole Blood",
    "Skin",
    "Brain / Muscle",
    "Esophagus / Skin",
    "Adipos / Blood Vessel / Esophagus / Skin",
    "Blood Vessel / Esophagus"
  )

gplots <- list()
for(m in 1:length(show_name)){

  ind = which(show_name[m] == terms_name)
  egene = terms_egene[ind,]
  name = terms_name[ind]
  pval = terms_pval[ind]
  pd = terms_pd[ind,]

  # color based on the relative scale of the eigengene
  # show the top 25% expressed gene
  colbar_len = 20
  colbar = colorRampPalette(c("black", "black", "black", "white"))(colbar_len)
  breaks = unique(quantile(egene, probs = seq(0, 1, length.out = colbar_len)))
  cut = cut(egene, breaks = breaks, include.lowest = TRUE)
  module_color = colbar[cut]

  df_emb <- data.frame(
    x = emb[,1],
    y = emb[,2]
  )
  scale_para <- (max(emb) - min(emb)) / 3
  df_arrow <- data.frame(
    xstart = 0,
    ystart = 0,
    xend = pd[1] * scale_para,
    yend = pd[2] * scale_para,
    label = name
  )

  # cur_name
  width = 30
  cur_name <- paste0(
    str_wrap(name, width = width),
    ifelse(nchar(name) > width, "\n(", "\n\n("),
    "P-value = ",
    sprintf("%.2e", pval),
    ")"
  )
  xlab <- paste0(
    str_wrap(show_label[m], width = width),
    ifelse(nchar(show_label[m]) > width, "", "\n")
  )

  g <- ggplot() +
    geom_point(data = df_emb,
               aes(x = x, y = y),
               color = module_color,
               size = 0.1) +
    geom_segment(data = df_arrow,
                 mapping = aes(x = xstart,
                               y = ystart,
                               xend = xend,
                               yend = yend),
                 col = "white",
                 arrow = arrow(length = unit(0.05, "npc")),
                 show.legend = FALSE) +
    ggtitle(cur_name) +
    labs(x = xlab, y = " ") +
    dark_theme_gray() +
    theme(plot.title = element_text(size = 14, family = "Arial"),
          axis.title = element_text(size = 14),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank())

  gplots[[m]] <- g
}

gridExtra::grid.arrange(grobs = gplots[1:4], nrow=1)
gridExtra::grid.arrange(grobs = gplots[5:8], nrow=1)


