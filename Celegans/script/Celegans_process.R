library(monocle3)
library(ggplot2)

## download data
expr_url <- url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds")
cell_meta_url <- url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds")
gene_meta_url <- url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds")

expr_mat <- readRDS(expr_url)
cell_meta <- readRDS(cell_meta_url)
gene_meta <- readRDS(gene_meta_url)

# filter 1: remove the genes with low expression
f1 <- rowSums(expr_mat) > 10
select_f1 <- names(f1)[f1]
expr_mat <- expr_mat[select_f1,]
gene_meta <- gene_meta[select_f1,]

# filter 2: remove the genes with low variance
rowvar <- rowVars(expr_mat)
f2 <- rowvar > quantile(rowvar, 0.25)
select_f2 <- names(f2)[f2]
expr_mat <- expr_mat[select_f2,]
gene_meta <- gene_meta[select_f2,]

# Batch correction (following Monocle 3)
cds <- new_cell_data_set(expr_mat,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_meta)
cds <- preprocess_cds(cds, num_dim = 50, method = "PCA")
cds <- align_cds(cds,
                 alignment_group = "batch",
                 residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
preprocess_mat <- SingleCellExperiment::reducedDims(cds)[["Aligned"]]

## Dimension Reduction with UMAP
umap_model <- uwot::umap(as.matrix(preprocess_mat),
                         n_components = 2,
                         metric = "cosine",
                         n_neighbors = 15,
                         min_dist = 0.5,
                         ret_model = TRUE)
set.seed(2024)
umap_res <- uwot::umap_transform(X=as.matrix(preprocess_mat), model=umap_model, n_threads=1)
emb <- scale(umap_res, center = TRUE, scale = FALSE)

exp <- as.matrix(t(cds@assays@data$counts))
exp <- exp[,match(colnames(exp), gene_meta$id)]
cell_meta <- cds@colData
colnames(exp) <- gene_meta$gene_short_name
cell_type <- cell_meta$cell.type
cell_type[is.na(cell_type)] <- "Unknown"
cell_meta$cell.type <- cell_type

saveRDS(emb, file = "./Celegans/results/global_emb.rds")
saveRDS(exp, file = "./Celegans/results/global_exp.rds")
saveRDS(cell_meta, file = "./Celegans/results/global_meta.rds")

## Download functional gene sets
out = "./Celegans/data/GO_Biological_Process_2018.txt"
url = "https://maayanlab.cloud/WormEnrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018"
download.file(url, destfile = out)







