library(ggplot2)
library(dplyr)

exp <- readRDS("./GTEx/results/GTEx_expr.rds")
emb    <- readRDS("./GTEx/results/GTEx_emb.rds")
meta   <- readRDS("./GTEx/results/GTEx_meta.rds")
tcol   <- readRDS("./GTEx/results/GTEx_tcol.rds")
mapped <- readRDS("./GTEx/results/mapped_genes.rds")

# focus on the two GO terms
gene_sets <- fgsea::gmtPathways("./GTEx/data/GO_Biological_Process_2023.txt")
term_name <- c("Striated Muscle Contraction", "Signal Release From Synapse")
term_go <- c("GO:0006941", "GO:0099643")

ind1 <- which(grepl(term_go[1], names(gene_sets)))
ind2 <- which(grepl(term_go[2], names(gene_sets)))
gsets <- gene_sets[c(ind1, ind2)]


# focus on the five genes
ShowGenes <- function (gsets, num) {
  gname <- colnames(exp)
  overlap_g  = intersect(gsets, gname)
  exp_ov     = exp[ , overlap_g]
  svd_res    = svd(exp_ov)
  cur_egene  = svd_res$u[,1]
  cur_egene  = as.numeric(ifelse(cor(cur_egene, rowMeans(exp_ov)) > 0, 1, -1))*cur_egene
  cur_pd     = ccaPP::ccaGrid(x = cur_egene, y = emb, method = "pearson")$B
  ng = length(overlap_g)
  coord = matrix(0, nrow = ng, ncol = 2)
  for(i in 1:ng){
    cca_res   = ccaPP::ccaGrid(exp_ov[,i], y = emb, method = "pearson")
    coord[i,] = cca_res$B
  }
  cossim = as.vector(coord %*% cur_pd)
  cossim = setNames(cossim, overlap_g)
  select_g <- overlap_g[order(cossim, decreasing = TRUE)[1:5]]

  gplots <- list()
  for(m in 1:num){

    gname = select_g[m]
    gexp = exp_ov[,gname]

    colbar_len = 20
    colbar = colorRampPalette(c("black", "black", "black", "white"))(colbar_len)
    breaks = unique(quantile(gexp, probs = seq(0, 1, length.out = colbar_len)))
    cut = cut(gexp, breaks = breaks, include.lowest = TRUE)
    module_color = colbar[cut]

    df_emb <- data.frame(
      x = emb[,1],
      y = emb[,2]
    )
    scale_para <- (max(emb) - min(emb)) / 3

    g <- ggplot() +
      geom_point(data = df_emb,
                 aes(x = x, y = y),
                 color = module_color,
                 size = 0.1) +
      ggtitle(gname) +
      labs(x = "", y = " ") +
      dark_theme_gray() +
      theme(plot.title = element_text(size = 14, family = "Arial"),
            axis.title = element_text(size = 14),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank())

    gplots[[m]] <- g
  }
  gplots
}

gplots1 <- ShowGenes(gsets$`Striated Muscle Contraction (GO:0006941)`, 4)
gridExtra::grid.arrange(grobs = gplots1, nrow=1)

gplots2 <- ShowGenes(gsets$`Signal Release From Synapse (GO:0099643)`, 4)
gridExtra::grid.arrange(grobs = gplots2, nrow=1)
