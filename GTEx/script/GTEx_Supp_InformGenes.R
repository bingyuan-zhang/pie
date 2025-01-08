source("./Rscript/functions.R")

## Input
emb <- readRDS(file="./GTEx/results/GTEx_emb.rds") # embedding matrix n*2
exp <- readRDS(file="./GTEx/results/GTEx_expr.rds") # expression matrix n*p
gene_sets <- fgsea::gmtPathways("./GTEx/data/GO_Biological_Process_2023.txt")

## Genes selected by PIE
# inputs
ns  = nrow(emb)
nd  = ncol(emb)
# parameters
method = "pearson"
ratio = 0.25
# Filter 1: SCmarker selection genes
rownames(exp) = 1:ns
obj_scm = SCMarker::ModalFilter(data=t(exp), geneK=10, cellK=10)
obj_scm = SCMarker::GeneFilter(obj=obj_scm)
obj_scm = SCMarker::getMarker(obj=obj_scm, n=10)
exp_scm = exp[, obj_scm$marker]
ng      = ncol(exp_scm)
# Map the selected genes
cor   = rep(0, ng)
coord = matrix(0, nrow = ng, ncol = nd)
for(i in 1:ng){
  cca_res   = ccaPP::ccaGrid(exp_scm[,i], y = emb, method = method)
  coord[i,] = cca_res$B
  cor[i]    = cca_res$cor
}
# select the top 75% correlated features regrading the embedding
topid  = which(cor > quantile(cor, ratio, na.rm = TRUE))
select_pie = obj_scm$marker[topid]

## Genes selected by SCMarker
select_scmarker <- obj_scm$marker

## HEG
num <- length(select_pie)
total_expression <- colSums(exp)
select <- order(total_expression, decreasing = TRUE)[1:num]
select_heg <- colnames(exp)[select]

## HVG
var_expression <- matrixStats::colVars(exp)
select <- order(var_expression, decreasing = TRUE)[1:num]
select_hvg <- colnames(exp)[select]

# show Venn Plot
library(VennDiagram)
areaA <- length(select_scmarker)
areaB <- length(select_pie)
areaC <- length(select_heg)
areaD <- length(select_hvg)

nAB <- length(intersect(select_scmarker, select_pie))
nAC <- length(intersect(select_scmarker, select_heg))
nAD <- length(intersect(select_scmarker, select_hvg))
nBC <- length(intersect(select_pie, select_heg))
nBD <- length(intersect(select_pie, select_hvg))
nCD <- length(intersect(select_heg, select_hvg))

nABC <- length(
  intersect(
    intersect(select_scmarker, select_pie),
    select_heg
  )
)
nABD <- length(
  intersect(
    intersect(select_scmarker, select_pie),
    select_hvg
  )
)
nACD <- length(
  intersect(
    intersect(select_scmarker, select_heg),
    select_hvg
  )
)
nBCD <- length(
  intersect(
    intersect(select_pie, select_heg),
    select_hvg
  )
)

nABCD <- length(
  Reduce(intersect, list(select_scmarker, select_pie, select_heg, select_hvg))
)

venn.plot <- draw.quad.venn(
  area1 = areaA,
  area2 = areaB,
  area3 = areaC,
  area4 = areaD,
  n12   = nAB,  n13   = nAC,  n14   = nAD,
  n23   = nBC,  n24   = nBD,  n34   = nCD,
  n123  = nABC, n124  = nABD, n134  = nACD, n234  = nBCD,
  n1234 = nABCD,


  category = c("SCMarker", "PIE", "HEG", "HVG"),
  fill = c("#99CCFF", "#FF9999", "#99FF99", "#FFFF99"),
  alpha = 0.6,
  cat.col = c("blue", "red", "darkgreen", "orange"),
  cat.cex = 1.2,
  cex = 1.2,
  main = "Venn Diagram"
)

grid.draw(venn.plot)
