library(ggplot2)
library(ggrepel)
library(ggdark)
library(purrr)
library(stringr)
library(dplyr)
library(matrixStats)

source("./Rscript/functions.R")
emb    <- readRDS("./GTEx/results/GTEx_emb_umap.rds")
meta   <- readRDS("./GTEx/results/GTEx_meta.rds")
tcol   <- readRDS("./GTEx/results/GTEx_tcol.rds")
mapped <- readRDS("./GTEx/results/mapped_genes.rds")
terms  <- readRDS("./GTEx/results/enriched_terms.rds")

## Enriched GO terms
terms_pd    <- terms$pd # projected direction (coordinates) of the eigengene
terms_padj  <- terms$padj # adjusted p-value of the GSEA test
terms_egene <- terms$egene # eigengene vector with length n
terms_name  <- terms$name %>% # name of the GO term
  strsplit(split = "\\(") %>%
  map_vec(~ head(.x, 1)) %>%
  str_remove("\\)") %>%
  str_remove(" $")

## Change the 2D coordinates of mapped GO terms into the polar coordinates
# Number of intervals
choose_k_interval <- function(k){
  nterm <- nrow(terms_pd)
  theta <- rep(0, nterm)
  for(i in 1:nterm) theta[i] <- coord2theta(terms_pd[i,])
  intv <- cut(theta, breaks=seq(0,2,length.out=k+1), include.lowest = TRUE)
  intv <- match(intv, levels(intv))
  
  ## Find the representative top GO terms for k intervals
  term_by_intv <- list()
  for(i in 1:k){
    
    ind = (intv==i)
    intv_name  = terms_name[ind]
    intv_padj  = terms_padj[ind]
    intv_pd    = terms_pd[ind,]
    intv_egene = terms_egene[ind,]
    
    if (length(intv_name) > 1) {
      
      reorder = order(intv_padj, decreasing = FALSE)
      intv_name = intv_name[reorder]
      intv_padj = intv_padj[reorder]
      intv_pd = intv_pd[reorder,]
      intv_egene = intv_egene[reorder,]
      
      term_by_intv[[paste("intv", i)]] <-
        list(name  = intv_name,
             padj  = intv_padj,
             pd    = intv_pd,
             egene = intv_egene)
    }
  }
  
  ## Figure 2A
  # add sample coordinates
  df_samples <- data.frame(
    x = emb[,1],
    y = emb[,2],
    smts = factor(meta$SMTS),
    smtsd = factor(meta$SMTSD)
  )
  lab <- matrix(0, nrow = 27, ncol=2)
  for(i in 1:27){
    idx <- meta$SMTSD == names(table(meta$SMTSD))[i]
    lab[i,] <- colMedians(emb[idx,])
  }
  # add sample labels
  df_sample_labels <- data_frame(
    x = lab[,1],
    y = lab[,2],
    label = names(table(meta$SMTSD))
  )
  
  ggplot_emb_samples <-
    ggplot() +
    geom_point(data = df_samples,
               aes(x=x, y=y),
               color=tcol$tcolor[factor(meta$SMTS)],
               alpha = 5/10,
               size = 0.05) +
    labs( x = "", y = "") +
    dark_theme_gray() +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text  = element_blank())
  
  
  ## Figure S2A
  # add mapped functions' coordinates
  scale_para <- (max(df_samples$x) - min(df_samples$x)) / 2
  df_func_coord <-
    data.frame(
      x = terms_pd[, 1] * scale_para,
      y = terms_pd[, 2] * scale_para,
      size = -log10(terms_padj) / max(-log10(terms_padj)) * 10
    )
  
  # add arrows
  select_num <- 3
  arrow_by_intv <- matrix(0, nrow=length(term_by_intv), ncol=2)
  for(i in 1:length(term_by_intv)) {
    x <- term_by_intv[[i]]
    if(length(x$padj) == 1) {
      arrow_by_intv[i,] <- x$pd
    } else {
      choose = min(select_num, nrow(x$pd))
      arrow_by_intv[i,] <- colMeans(x$pd[1:choose, ])
    }
  }
  df_func_arrow <-
    data.frame(
      xstart = rep(0, length(term_by_intv)),
      ystart = rep(0, length(term_by_intv)),
      xend = arrow_by_intv[,1]*scale_para * 0.9,
      yend = arrow_by_intv[,2]*scale_para * 0.9
    )
  
  # add terms labels
  df_func_label <- data.frame()
  for(i in 1:length(term_by_intv)){
    ind <- 1:min(select_num,nrow(term_by_intv[[i]]))
    label <- term_by_intv[[i]]$name[ind]
    x <- term_by_intv[[i]]$pd[,1][ind]*1.5*scale_para
    y <- term_by_intv[[i]]$pd[,2][ind]*1.5*scale_para
    df_func_label <- rbind(df_func_label, data.frame(x=x, y=y, label=label))
  }
  
  set.seed(2024)
  figure_c <-
    ggplot_emb_samples +
    geom_point(data = df_func_coord,
               aes(x = x, y = y),
               color = "white",
               size = df_func_coord$size,
               show.legend = FALSE) +
    geom_segment(data = df_func_arrow,
                 mapping = aes(x = xstart,
                               y = ystart,
                               xend = xend,
                               yend = yend),
                 arrow = arrow(length = unit(0.05, "npc")),
                 size = 1.5,
                 color = "white",
                 alpha = 7/10,
                 show.legend = FALSE) +
    geom_label_repel(data = df_func_label,
                     aes(label=str_wrap(label,20), x=x, y=y),
                     family = "Arial",
                     box.padding = 0.1,
                     point.padding = 0.1,
                     segment.size = 0,
                     force = 1,
                     force_pull = 0.5,
                     size = 9/2.845,
                     colour="white",
                     fill = alpha("white",0),
                     max.overlaps = 30)
  figure_c
}


figure_S2A <- choose_k_interval(4)
figure_S2B <- choose_k_interval(8)
figure_S2C <- choose_k_interval(24)

gridExtra::grid.arrange(grobs = list(figure_S2A, figure_S2B, figure_S2C), nrow=1)
