library(ComplexHeatmap)
library(dplyr)
library(corrplot)

emb    <- readRDS("./GTEx/results/GTEx_emb.rds")
meta   <- readRDS("./GTEx/results/GTEx_meta.rds")
tcol   <- readRDS("./GTEx/results/GTEx_tcol.rds")
terms  <- readRDS("./GTEx/results/enriched_terms.rds")
term_by_intv <- readRDS("./GTEx/results//term_intv.rds")

## Heatmap
ht_opt$message = TRUE

# prepare eigenegene matrix for functions
select_num <- 3
egenesmat <- c()
for(i in 1:length(term_by_intv)){
  mat_add <- term_by_intv[[i]]$egene[1:select_num, ]
  rownames(mat_add) <- term_by_intv[[i]]$name[1:select_num]
  egenesmat <- rbind(egenesmat, mat_add)
}

# column annotation
anno_column <- HeatmapAnnotation(
    ts = meta$SMTS,
    ts_sub = meta$SMTSD,
    col = list(ts = tcol$tcolor, ts_sub = tcol$tscolor),
    annotation_legend_param = list(
      ts = list(
        title = "Tissue",
        grid_height = unit(6, "mm"),
        grid_width = unit(6, "mm") #,
        # title_gp = gpar(fontsize = 14, fontface = "bold"),
        # labels_gp = gpar(fontsize = 12)
      ),
      ts_sub = list(
        title = "Subcategory",
        grid_height = unit(6, "mm"),
        grid_width = unit(6, "mm") #,
        # title_gp = gpar(fontsize = 14, fontface = "bold"),
        # labels_gp = gpar(fontsize = 12)
      )
    ),
    show_annotation_name = c(tissue = FALSE, sub = FALSE),
    show_legend = c(tissue = FALSE, sub = FALSE)
  )

# heatmap parameters
heatmap_params <- list(
  mat = egenesmat,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(20),
  breaks = unique(quantile(as.matrix(egenesmat), probs = seq(0,1,length.out=20))),
  annotation_col = anno_column
  )

# make heatmap object
figure_heatmap <-
  ComplexHeatmap::pheatmap(
    mat = heatmap_params$mat,
    color = heatmap_params$color,
    breaks = heatmap_params$breaks,
    bottom_annotation = heatmap_params$annotation_col,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    row_names_max_width = unit(20, "cm"),
    fontsize = 11,
    fontfamily = "ArialMT",
    use_raster = FALSE,
    legend = FALSE
  )

draw(figure_heatmap, padding = unit(c(2, 2, 2, 2), "cm"))

################################################################################
Mcor <- cor(t(egenesmat))
par(family = "ArialMT")
corrplot(
  Mcor,
  method = 'square',
  order = 'hclust',
  cl.pos = 'b',
  tl.pos = 'l',
  tl.col = 'black',
  tl.cex = 1
)




