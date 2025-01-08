library(SCMarker)
library(fgsea)

MapGenes <- function(inputs) {
  # inputs
  emb = inputs$emb
  exp = inputs$exp
  ns  = nrow(emb)
  nd  = ncol(emb)

  # parameters
  method = "pearson"
  ratio = 0.25

  # Filter 1: SCmarker selection genes
  print(paste("Filtering the genes"))
  tic <- proc.time()
  rownames(exp) = 1:ns
  obj_scm = SCMarker::ModalFilter(data=t(exp), geneK=10, cellK=10)
  obj_scm = SCMarker::GeneFilter(obj=obj_scm)
  obj_scm = SCMarker::getMarker(obj=obj_scm, n=10)
  exp_scm = exp[, obj_scm$marker]
  ng      = ncol(exp_scm)
  toc <- proc.time()
  time_scmarker <- (toc - tic)[3]

  # Map the selected genes
  print(paste("Mapping the genes"))
  tic <- proc.time()
  cor   = rep(0, ng)
  coord = matrix(0, nrow = ng, ncol = nd)
  for(i in 1:ng){
    cca_res   = ccaPP::ccaGrid(exp_scm[,i], y = emb, method = method)
    coord[i,] = cca_res$B
    cor[i]    = cca_res$cor
  }
  toc <- proc.time()
  time_map <- (toc - tic)[3]

  # select the top 75% correlated features regrading the embedding
  topid  = which(cor > quantile(cor, ratio, na.rm = TRUE))
  coord  = coord[topid, ]
  cor    = cor[topid]
  select = obj_scm$marker[topid]
  rownames(coord) = names(cor) = select

  res = list(select = select,
             coord = coord,
             cor = cor,
             time_scmarker = time_scmarker,
             time_map = time_map,
             inputs_enrich_test = list(exp = exp[,select],
                                       emb = emb,
                                       coord = coord)
             )
  return(res)
}

EnrichTest <- function(inputs, gene_sets){
  # inputs
  exp     = inputs$exp
  emb     = inputs$emb
  coord   = inputs$coord
  gname   = colnames(inputs$exp)
  gs_num  = length(gene_sets)
  gs_name = names(gene_sets)

  # parameters
  p_thres  = 0.05
  min_olap = 3
  max_olap = 300
  remind   = 500

  term_name  = c()
  term_pval  = c()
  term_nes  = c()
  term_pd    = c()
  term_egene = c()
  time_svd = 0
  time_map = 0
  time_gsea = 0

  for(i in 1:gs_num){
    if(i%%remind == 0) print(paste0("Testing the ", i, "-th GO term."))

    term_gs_name = gene_sets[[i]]
    overlap_g  = intersect(term_gs_name, gname)

    # check conditions (the overlapped gene set ranges)
    cond1 = length(overlap_g) <= min_olap
    cond2 = length(overlap_g) > max_olap
    if(cond1 | cond2) next

    # extract the eigengene
    tic <- proc.time()
    exp_ov     = exp[ , overlap_g]
    svd_res    = svd(exp_ov)
    cur_egene  = svd_res$u[,1]
    cur_egene  = as.numeric(ifelse(cor(cur_egene, rowMeans(exp_ov)) > 0, 1, -1))*cur_egene
    toc <- proc.time()
    time_svd <- time_svd + (toc - tic)[3]

    # map the eigengene
    tic <- proc.time()
    cur_pd     = ccaPP::ccaGrid(x = cur_egene, y = emb, method = "pearson")$B
    toc <- proc.time()
    time_map <- time_map + (toc - tic)[3]

    # rank the mapped genes
    tic <- proc.time()
    cossim = as.vector(coord %*% cur_pd)
    cossim = setNames(cossim, gname)
    grank  = rank(cossim, ties.method = "first")
    grank  = na.omit(grank)
    grank  = grank[unique(names(grank))]

    # perform GSEA test
    gsea_res = fgsea::fgseaMultilevel(pathways = gene_sets[i],
                                      stats = grank,
                                      minSize = min_olap,
                                      maxSize = max_olap,
                                      scoreType = "pos",
                                      gseaParam = 1)
    toc <- proc.time()
    time_gsea <- time_gsea + (toc - tic)[3]

    cur_pval = gsea_res$pval
    cur_nes = gsea_res$NES

    if(length(cur_pval) != 0) {
      term_name  = c(term_name, gs_name[i])
      term_pval  = c(term_pval, cur_pval)
      term_nes   = c(term_nes, cur_nes)
      term_pd    = cbind(term_pd, cur_pd)
      term_egene = cbind(term_egene, cur_egene)
    }
  }

  ind   = (term_pval < p_thres)
  name  = term_name[ind]
  pval  = term_pval[ind]
  nes   = term_nes[ind]
  pd    = t(term_pd)[ind,]
  egene = t(term_egene)[ind,]

  list(
    name = name,
    pd = pd,
    pval = pval,
    nes = nes,
    egene = egene,
    time_svd = time_svd,
    time_map = time_map,
    time_gsea = time_gsea
  )
}

coord2theta <- function(coord){
  x = coord[1]
  y = coord[2]

  if(x>=0 & y>=0) {
    asin(y)/pi
  }else if(x<0 & y>=0){
    1 - asin(y)/pi
  }else if(x<0 & y<0){
    1 + asin(-y)/pi
  }else {
    2 - asin(-y)/pi
  }
}

ggColorHue <- function(n, l=65) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=100)[1:n]
}


