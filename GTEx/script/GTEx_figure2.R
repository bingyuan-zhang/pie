library(ggplot2)
library(ggrepel)
library(ggdark)
library(purrr)
library(stringr)
library(dplyr)

source("./R/functions.R")
emb    <- readRDS("./GTEx/results/GTEx_emb.rds")
meta   <- readRDS("./GTEx/results/GTEx_meta.rds")
tcol   <- readRDS("./GTEx/results/GTEx_tcol.rds")
mapped <- readRDS("./GTEx/results/mapped_genes.rds")
terms  <- readRDS("./GTEx/results/enriched_terms.rds")

## Enriched GO terms
terms_pd    <- terms$pd # projected direction (coordinates) of the eigengene
terms_pval  <- terms$pval # p-value of the GSEA test
terms_egene <- terms$egene # eigengene vector with length n
terms_name  <- terms$name %>% # name of the GO term
  strsplit(split = "\\(") %>%
  map_vec(~ head(.x, 1)) %>%
  str_remove("\\)") %>%
  str_remove(" $")

## Change the 2D coordinates of mapped GO terms into the polar coordinates
## Divide 12 equally spaced intervals by angle
nterm <- nrow(terms_pd)
theta <- rep(0, nterm)
for(i in 1:nterm) theta[i] <- coord2theta(terms_pd[i,])
k <- 12
intv <- cut(theta, breaks=seq(0,2,length.out=k+1), include.lowest = TRUE)
intv <- match(intv, levels(intv))

## Find the representative top GO terms for 12 intervals
term_by_intv <- list()
for(i in 1:k){
  ind = (intv==i)
  intv_name  = terms_name[ind]
  intv_pval  = terms_pval[ind]
  intv_pd    = terms_pd[ind,]
  intv_egene = terms_egene[ind,]

  if (length(intv_name) != 0) {

    reorder = order(intv_pval, decreasing = FALSE)
    intv_name = intv_name[reorder]
    intv_pval = intv_pval[reorder]
    intv_pd = intv_pd[reorder,]
    intv_egene = intv_egene[reorder,]

    term_by_intv[[paste("intv", i)]] <-
      list(name  = intv_name,
           pval  = intv_pval,
           pd    = intv_pd,
           egene = intv_egene)
  }
}
saveRDS(term_by_intv, file="./GTEx/results//term_intv.rds")

## Figure 2A
term_by_intv <- readRDS("./GTEx/results//term_intv.rds")
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

set.seed(2024)
figure_a <-
  ggplot_emb_samples +
  geom_text_repel(data = df_sample_labels,
                  aes(x=x,y=y, label=str_wrap(label, 20)),
                  size = 13/2.845,
                  segment.size = 0.3,
                  box.padding = 0.2,
                  family = "Arial",
                  max.overlaps = 10)
figure_a

## Figure 2B
k = 12
select_num = 3
gene_num <- length(mapped$select)
theta_gene <- c()
for(i in 1:gene_num) theta_gene[i] <- coord2theta(mapped$coord[i,])
intv1 <- cut(theta_gene, breaks=seq(0,2,length.out=k+1), include.lowest = TRUE)
intv1 <- match(intv1, levels(intv1))
scale_para <- (max(df_samples$x) - min(df_samples$x)) / 2.5 / 5
df_gene_coord <- data.frame(
  x = mapped$coord[,1] * (mapped$cor + 4.5) * scale_para,
  y = mapped$coord[,2] * (mapped$cor + 4.5) * scale_para,
  col = factor(intv1)
)
df_gene_labels <- c()
for(i in 1:k){
  idx <- which(intv1==i)
  name <- mapped$select[idx]
  idxx <- order(mapped$cor[idx], decreasing = TRUE)[1:select_num]
  name <- name[idxx]
  idxx <- idx[idxx]
  df_gene_labels <- rbind(df_gene_labels,
                          data.frame(x = mapped$coord[idxx, 1],
                                     y = mapped$coord[idxx, 2],
                                     label = name))
}

set.seed(2024)
figure_b <-
  ggplot_emb_samples +
  geom_point(data = df_gene_coord,
             aes(x = x, y = y),
             size = 0.2,
             show.legend = FALSE)  +
  geom_label_repel(data = df_gene_labels,
                   aes(label=str_wrap(label), x=x*18, y=y*18),
                   box.padding = 0.1,
                   segment.size = 0,
                   force_pull = 0.1,
                   size = 13/2.845,
                   colour="white",
                   fill = alpha("white",0),
                   family = "Arial",
                   max.overlaps = 50)
figure_b


## Figure 2C
# add mapped functions' coordinates
scale_para <- (max(df_samples$x) - min(df_samples$x)) / 2
df_func_coord <-
  data.frame(
    x = terms_pd[, 1] * scale_para,
    y = terms_pd[, 2] * scale_para,
    size = -log10(terms_pval) / max(-log10(terms_pval)) * 10
  )

# add arrows
select_num <- 3
arrow_by_intv <- matrix(0, nrow=length(term_by_intv), ncol=2)
for(i in 1:length(term_by_intv)) {
  x <- term_by_intv[[i]]
  if(length(x$pval) == 1) {
    arrow_by_intv[i,] <- x$pd
  } else {
    choose = min(select_num, nrow(x$pd))
    arrow_by_intv[i,] <- colMeans(x$pd[1:choose, ])
  }
}
df_func_arrow <-
  data.frame(
    xstart = rep(0, 9),
    ystart = rep(0, 9),
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


save(figure_a, figure_b, figure_c, file="./GTEx/results/Figure2.rdata")
