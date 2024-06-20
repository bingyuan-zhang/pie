library(ggplot2)
library(ggrepel)
library(ggdark)
library(purrr)
library(stringr)
library(dplyr)

source("./R/functions.R")

## Find representative terms
# enriched GO terms
terms <- readRDS("./Celegans/results/global_enriched_terms.rds")
terms_pd <- terms$pd
terms_pval <- terms$pval
terms_egene <- terms$egene
terms_name <- gsub(" \\(GO:[0-9]+\\)", "", terms$name)
terms_id <- sub(".*\\(GO:([0-9]+)\\)", "GO:\\1", terms$name)

term_num <- nrow(terms_pd)
theta <- rep(0, length=term_num)
for(i in 1:term_num) theta[i] <- coord2theta(terms_pd[i,])
k = 12
intv <- cut(theta, breaks=seq(0,2,length.out=k+1), include.lowest = TRUE)
intv <- match(intv, levels(intv))
term_by_intv <- list()
for(i in 1:k){
  cur_name  <- terms_name[intv==i]
  cur_pval  <- terms_pval[intv==i]
  cur_pd    <- terms_pd[intv==i,]
  cur_egene <- terms_egene[intv==i,]
  if(length(cur_pval) == 0) {
    next
  } else if(length(cur_pval) > 1){
    idx <- order(cur_pval, decreasing = FALSE)
    term_by_intv[[i]] <-
      list(name = cur_name[idx],
           pval = cur_pval[idx],
           pd = cur_pd[idx,],
           egene = cur_egene[idx,])
  } else {
    term_by_intv[[i]] <-
      list(name = cur_name,
           pval = cur_pval,
           pd = cur_pd,
           egene = cur_egene)
  }
}

## figure 6A
emb <- readRDS("./Celegans/results/global_emb.rds")
meta <- readRDS("./Celegans/results/global_meta.rds")
cell_type <- meta$cell.type
embryo_time <- meta$embryo.time
ncolor <- ggColorHue(28)
names(ncolor) <- unique(cell_type)
df_emb <- data.frame(
  x = emb[,1],
  y = emb[,2],
  cell_type = cell_type,
  cell_col = ncolor[cell_type],
  time = embryo_time
)
cell_type_num <- length(table(cell_type))
lab <- matrix(0, nrow = cell_type_num, ncol=2)
for(i in 1:cell_type_num){
  idx <- names(table(cell_type))[i] == cell_type
  lab[i,] <- colMedians(emb[idx,])
}
ind <- nchar(names(table(cell_type))) <= 4
df_emb_text <- data_frame(
  x = lab[ind,1],
  y = lab[ind,2],
  label = names(table(cell_type))[ind]
)


figure_6A <-
  ggplot() +
  geom_point(data = df_emb,
             aes(x = x,
                 y = y,
                 color = time),
             show.legend = FALSE,
             size = 0.1) +
  labs(x = "", y = "") +
  geom_text_repel(data = df_emb_text,
                  family = "Arial",
                  aes(x = x,
                      y = y,
                      label = label),
                  size = 18/2.854,
                  max.overlaps = 50) +
  dark_theme_grey() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

figure_6A

## figure 6B
scale_para <- (max(df_emb$x) - min(df_emb$x)) / 2
df_mapped_func <- data.frame(
  x = terms_pd[,1]*scale_para,
  y = terms_pd[,2]*scale_para,
  size = -log10(terms_pval)/ max(-log10(terms_pval)) * 7
)
# mapped arrows
select_num <- 3
Arrows <- c()
for(i in 1:k){
  # print(i)
  pds <- term_by_intv[[i]]$pd
  num <- min(select_num, nrow(pds))
  if(length(pds) < 2) {
    next
  } else if(length(pds) == 2) {
    Arrows <- rbind(Arrows, pds)
  } else {
    Arrows <- rbind(Arrows, colMeans(pds[1:num, ]))
  }
}
df_arrow_func <- data.frame(
  xstart = rep(0, nrow(Arrows)),
  ystart = rep(0, nrow(Arrows)),
  xend = Arrows[,1]*scale_para * 0.9,
  yend = Arrows[,2]*scale_para * 0.9,
  yend = Arrows[,2]*scale_para * 0.9
)
# mapped labels
df_lab <- data.frame()
j = 0
for(i in 1:k){
  if(is.null(term_by_intv[[i]])) next
  j = j + 1
  num <- min(select_num, length(term_by_intv[[i]]$pval))
  label <- term_by_intv[[i]]$name[1:num]
  xend <- rep(Arrows[j,1]*1.7*scale_para, num)
  yend <- rep(Arrows[j,2]*1.7*scale_para, num)
  col <- rep(color[i], num)

  df_lab <- rbind(df_lab, data.frame(x=xend, y=yend, label=label, col=col))
}

set.seed(2024)
figure_6B <- ggplot() +
  geom_point(data = df_emb,
             aes(x = x,
                 y = y),
             col = df_emb$cell_col,
             alpha = 0.4,
             show.legend = FALSE,
             size = 0.1) +
  geom_point(data = df_mapped_func,
             aes(x=x, y=y),
             color = "white",
             size = df_mapped_func$size) +
  geom_segment(data = df_arrow_func,
               mapping = aes(x = xstart,
                             y = ystart,
                             xend = xend,
                             yend = yend),
               arrow = arrow(length = unit(0.05, "npc")),
               size = 1.5,
               color = alpha("white", 0.8),
               show.legend = FALSE) +
  geom_label_repel(data = df_lab,
                   aes(label=str_wrap(label, 18),
                       x = x,
                       y = y),
                   family = "Arial",
                   box.padding = 0.1,
                   segment.size = 0,
                   force_pull = 0.13,
                   force = 0.15,
                   size = 13/2.854,
                   color="white",
                   fill = alpha("white",0),
                   max.overlaps = 40) +
  dark_theme_gray() +
  labs(x = "", y = "") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

figure_6B

## figure 6C
exp <- readRDS("./Celegans/results/global_exp.rds")
mapped <- readRDS("./Celegans/results/global_mapped_genes.rds")
gsets <- fgsea::gmtPathways("./Celegans/data/GO_Biological_Process_2018.txt")

show.term <- function(ind) {
  name <- terms_name[ind]
  pval <- terms_pval[ind]
  id_gmt <- grep(terms_id[ind], names(gsets))
  idx <- which(colnames(exp) %in% gsets[[id_gmt]] & colnames(exp) %in% mapped$select)
  ovlp_g <- colnames(exp)[idx]

  gexpr_module <- exp[,idx]
  svd_res <- svd(gexpr_module)
  egene <- svd_res$u[,1]
  egene <- as.numeric(ifelse(cor(egene, rowMeans(gexpr_module)) > 0, 1, -1))*egene

  colbar_len = 20
  colbar <- colorRampPalette(c("black", rep("white", 4) ))(colbar_len)
  breaks <- seq(min(egene), max(egene), length.out = colbar_len)
  cut <- cut(egene, breaks = breaks, include.lowest = TRUE)
  module_color <- colbar[cut]

  df <- data.frame(
    x = emb[,1],
    y = emb[,2]
  )
  scale_para <- (max(emb) - min(emb)) / 3
  m <- which(terms_name == name)
  df_arrow <- data.frame(
    xstart = 0,
    ystart = 0,
    xend = terms_pd[m, 1] * scale_para,
    yend = terms_pd[m, 2] * scale_para,
    label = name
  )

  width = 30
  cur_name <- paste0(
    str_wrap(name, width = width),
    ifelse(nchar(name) > width, "\n(", "\n\n("),
    "P-value = ",
    sprintf("%.2e", pval),
    ")"
  )

  g <- ggplot() +
    geom_point(data = df,
               aes(x = x,
                   y = y),
               colour = alpha(module_color, 0.7),
               size = 0.1) +
    geom_segment(data = df_arrow,
                 mapping = aes(x = xstart,
                               y = ystart,
                               xend = xend,
                               yend = yend),
                 color = "white",
                 arrow = arrow(length = unit(0.05, "npc")),
                 show.legend = FALSE,
                 size = 1.5) +
    ggtitle(cur_name) +
    labs(x = "AFD branch", y = "") +
    dark_theme_gray() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 20),
          title = element_text(size = 18))
  g
  res <- list(g = g, overlap_gname = ovlp_g)
  return(res)
}

ind <- which(terms_name == "cGMP metabolic process")
res <- show.term(ind)
figure_6C <- res$g
figure_6C

## figure 6D
show.gene <- function(name){
  idx <- which(name == colnames(exp))
  gene <- exp[,idx]

  colbar_len = 20
  colbar <- colorRampPalette(c("black", rep("white", 4) ))(colbar_len)
  breaks <- seq(min(gene), max(gene), length.out = colbar_len)
  cut <- cut(gene, breaks = breaks, include.lowest = TRUE)
  module_color <- colbar[cut]

  df_emb <- data.frame(x = emb[,1], y = emb[,2])
  scale.para <- (max(emb) - min(emb))/5
  g <- ggplot() +
    geom_point(data = df_emb,
               aes(x = x,
                   y = y),
               color = alpha(module_color, 0.7),
               size = 0.1) +
    ggtitle(name) +
    labs(x = "", y = "") +
    dark_theme_gray() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          title = element_text(size = 17))
  g
}

ind <- which(terms_name == "cGMP metabolic process")
res <- show.term(ind)
glist <- list()
for(i in 1:length(res$overlap_gname)){
  glist[[i]] <- show.gene(res$overlap_gname[i])
}
figure_6D <- gridExtra::grid.arrange(grobs = glist, nrow=2)
figure_6D




