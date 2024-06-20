library(ggplot2)
library(ggrepel)
library(ggdark)
library(purrr)
library(stringr)

source("./R/functions.R")
gsets <- fgsea::gmtPathways("./Celegans/data/GO_Biological_Process_2018.txt")

## enriched global GO terms
terms_local <- readRDS("./Celegans/results/local_enriched_terms.rds")
terms_local_pd <- terms_local$pd
terms_local_pval <- terms_local$pval
terms_local_egene <- terms_local$egene
terms_local_name <- gsub(" \\(GO:[0-9]+\\)", "", terms_local$name)
terms_local_id <- sub(".*\\(GO:([0-9]+)\\)", "GO:\\1", terms_local$name)

emb_local <- readRDS("./Celegans/results/local_emb.rds")
exp_local <- readRDS("./Celegans/results/local_exp.rds")
mapped_local <- readRDS("./Celegans/results/local_mapped_genes.rds")

## enriched local GO terms
terms_global <- readRDS("./Celegans/results/global_enriched_terms.rds")
terms_global_pd <- terms_global$pd
terms_global_pval <- terms_global$pval
terms_global_egene <- terms_global$egene
terms_global_name <- gsub(" \\(GO:[0-9]+\\)", "", terms_global$name)
terms_global_id <- sub(".*\\(GO:([0-9]+)\\)", "GO:\\1", terms_global$name)

emb_global <- readRDS("./Celegans/results/global_emb.rds")
exp_global <- readRDS("./Celegans/results/global_exp.rds")
mapped_global <- readRDS("./Celegans/results/global_mapped_genes.rds")

## find representative local enriched GO terms
term_local_num <- nrow(terms_local_pd)
theta <- rep(0, length = term_local_num)
for(i in 1:term_local_num) theta[i] <- coord2theta(terms_local_pd[i,])
k = 12
intv <- cut(theta, breaks=seq(0,2,length.out=k+1), include.lowest = TRUE)
intv <- match(intv, levels(intv))
term_by_intv_local <- list()
for(i in 1:k){
  cur_name  <- terms_local_name[intv==i]
  cur_pval  <- terms_local_pval[intv==i]
  cur_pd    <- terms_local_pd[intv==i,]
  cur_egene <- terms_local_egene[intv==i,]
  if(length(cur_pval) == 0) {
    next
  } else if(length(cur_pval) > 1){
    idx <- order(cur_pval, decreasing = FALSE)
    term_by_intv_local[[i]] <-
      list(name = cur_name[idx],
           pval = cur_pval[idx],
           pd = cur_pd[idx,],
           egene = cur_egene[idx,])
  } else {
    term_by_intv_local[[i]] <-
      list(name = cur_name,
           pval = cur_pval,
           pd = cur_pd,
           egene = cur_egene)
  }
}

## figure 7A
emb_local <- readRDS("./Celegans/results/local_emb.rds")
meta_local <- readRDS("./Celegans/results/local_meta.rds")
cell_type_local <- meta_local$cell.type
embryo_time_local <- meta_local$embryo.time

meta_global <- readRDS("./Celegans/results/global_meta.rds")
cell_type_global <- meta_global$cell.type
ncolor <- ggColorHue(28)
names(ncolor) <- unique(cell_type_global)

df_emb_local <- data.frame(
  x = emb_local[,1],
  y = emb_local[,2],
  cell_type = cell_type_local,
  cell_col = ncolor[cell_type_local],
  time = embryo_time_local
)

cell_type_local_num <- length(table(cell_type_local))
lab <- matrix(0, nrow = cell_type_local_num, ncol=2)
for(i in 1:cell_type_local_num){
  idx <- names(table(cell_type_local))[i] == cell_type_local
  lab[i,] <- colMedians(emb_local[idx,])
}
df_emb_local_text <- data_frame(
  x = lab[,1],
  y = lab[,2],
  label = names(table(cell_type_local))
)

scale_para <- (max(df_emb_local$x) - min(df_emb_local$x)) / 1.5
df_mapped_func_local <- data.frame(
  x = terms_local_pd[,1]*scale_para,
  y = terms_local_pd[,2]*scale_para,
  size = -log10(terms_local_pval)/ max(-log10(terms_local_pval)) * 7
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
df_arrow_func_local <- data.frame(
  xstart = rep(0, nrow(Arrows)),
  ystart = rep(0, nrow(Arrows)),
  xend = Arrows[,1]*scale_para * 0.9,
  yend = Arrows[,2]*scale_para * 0.9,
  yend = Arrows[,2]*scale_para * 0.9
)

# mapped labels
df_lab_local <- data.frame()
j = 0
for(i in 1:k){
  if(is.null(term_by_intv_local[[i]])) next
  j = j + 1
  num <- min(select_num, length(term_by_intv_local[[i]]$pval))
  label <- term_by_intv_local[[i]]$name[1:num]
  xend <- rep(Arrows[j,1]*1.5*scale_para, num)
  yend <- rep(Arrows[j,2]*1.5*scale_para, num)
  df_lab_local <- rbind(df_lab_local,
                        data.frame(x=xend, y=yend, label=label))
}

set.seed(2024)
figure_7A <- ggplot() +
  geom_point(data = df_emb_local,
             aes(x = x,
                 y = y),
             col = df_emb_local$cell_col,
             alpha = 0.6,
             show.legend = FALSE,
             size = 0.7) +
  geom_text_repel(data = df_emb_local_text,
                  family = "Arial",
                  aes(x = x, y = y, label = label),
                  size = 15/2.854,
                  max.overlaps = 50) +
  geom_point(data = df_mapped_func_local,
             aes(x=x, y=y),
             color = "white",
             size = df_mapped_func_local$size) +
  geom_segment(data = df_arrow_func_local,
               mapping = aes(x = xstart,
                             y = ystart,
                             xend = xend,
                             yend = yend),
               arrow = arrow(length = unit(0.05, "npc")),
               size = 1.5,
               color = alpha("white", 0.6)) +
  geom_label_repel(data = df_lab_local,
                   aes(label=str_wrap(label, 20),
                       x = x,
                       y = y),
                   family = "Arial",
                   box.padding = 0.1,
                   segment.size = 0,
                   force_pull = 0.1,
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

figure_7A

## figure 7B
show.term.global <- function(ind) {
  name <- terms_global_name[ind]
  pval <- terms_global_pval[ind]
  id_gmt <- grep(terms_global_id[ind], names(gsets))
  idx <- which(colnames(exp_global) %in% gsets[[id_gmt]] &
                 colnames(exp_global) %in% mapped_global$select)
  ovlp_g <- colnames(exp_global)[idx]
  gexpr <- exp_global[,idx]

  svd_res <- svd(gexpr)
  egene <- svd_res$u[,1]
  egene <- as.numeric(ifelse(cor(egene, rowMeans(gexpr)) > 0, 1, -1))*egene

  colbar_len = 20
  colbar <- colorRampPalette(c("black", rep("white", 4) ))(colbar_len)
  breaks <- seq(min(egene), max(egene), length.out = colbar_len)
  cut <- cut(egene, breaks = breaks, include.lowest = TRUE)
  module_color <- colbar[cut]

  df <- data.frame(
    x = emb_global[,1],
    y = emb_global[,2]
  )
  scale_para <- (max(emb_global) - min(emb_global)) / 3
  m <- which(terms_global_name == name)
  df_arrow <- data.frame(
    xstart = 0,
    ystart = 0,
    xend = terms_global_pd[m, 1] * scale_para,
    yend = terms_global_pd[m, 2] * scale_para,
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
               colour = alpha(module_color, 0.5),
               size = 0.5) +
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
    labs(x = "Global View", y = "") +
    dark_theme_gray() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 20),
          title = element_text(size = 18))
  g
  res <- list(g = g, overlap_gname = ovlp_g)
  return(res)
}

show.term.local <- function(ind){
  name <- terms_local_name[ind]
  pval <- terms_local_pval[ind]
  id_gmt <- grep(terms_local_id[ind], names(gsets))
  idx <- which(colnames(exp_local) %in% gsets[[id_gmt]] &
                 colnames(exp_local) %in% mapped_local$select)
  ovlp_g <- colnames(exp_local)[idx]
  gexpr <- exp_local[,idx]

  svd_res <- svd(gexpr)
  egene <- svd_res$u[,1]
  egene <- as.numeric(ifelse(cor(egene, rowMeans(gexpr)) > 0, 1, -1))*egene

  colbar_len = 20
  colbar <- colorRampPalette(c("black", rep("white", 4) ))(colbar_len)
  breaks <- seq(min(egene), max(egene), length.out = colbar_len)
  cut <- cut(egene, breaks = breaks, include.lowest = TRUE)
  module_color <- colbar[cut]

  df <- data.frame(
    x = emb_local[,1],
    y = emb_local[,2]
  )
  scale_para <- (max(emb_local) - min(emb_local)) / 3
  m <- which(terms_local_name == name)
  df_arrow <- data.frame(
    xstart = 0,
    ystart = 0,
    xend = terms_local_pd[m, 1] * scale_para,
    yend = terms_local_pd[m, 2] * scale_para,
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
               colour = alpha(module_color, 0.5),
               size = 0.5) +
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
    labs(x = "Local View", y = "") +
    dark_theme_gray() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 20),
          title = element_text(size = 18))
  g
  res <- list(g = g, overlap_gname = ovlp_g)
  return(res)
}

ind_local <- which(terms_local_name == "cilium assembly")
g1 <- show.term.local(ind_local)$g
ind_global <- which(terms_global_name == "cilium assembly")
g2 <- show.term.global(ind_global)$g
ind_local <- which(terms_local_name == "mRNA splicing, via spliceosome")
g3 <- show.term.local(ind_local)$g
ind_global <- which(terms_global_name == "mRNA splicing, via spliceosome")
g4 <- show.term.global(ind_global)$g
gridExtra::grid.arrange(grobs = list(g1, g2, g3, g4), nrow=1)

