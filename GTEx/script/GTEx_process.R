# Load Packages
library(data.table)
library(dplyr)
library(parallel)
source("./R/functions.R")

# Load GTEx dataset and Meta data
url_gtpm <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
url_meta_sample <- "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
gtpm <- fread(url_gtpm, data.table = FALSE)
meta <- fread(url_meta_sample)

# out = "./GTEx/data/protein-coding_gene.txt"
# url = "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt"
# download.file(url, destfile = out)
gprotein <- fread("./GTEx/data/protein-coding_gene.txt", data.table = FALSE)

# focus on the protein coding genes
gtpm <- gtpm %>% filter(Description %in% gprotein$symbol)

# remove duplicated gene symbols
distinct_id <- match(unique(gtpm$Description), gtpm$Description)
gtpm <- gtpm[distinct_id, ]

# gene name as rownames
gname <- gtpm %>% pull(Description)
rownames(gtpm) <- gname

# map the meta data
gtpm <- gtpm[, colnames(gtpm) %in% meta$SAMPID]
idx <- match(colnames(gtpm), meta$SAMPID)
meta <- meta[idx, ]

# keep the high quality samples (SMRIN > 7)
hqc_idx <- meta$SMRIN > 7
gtpm <- gtpm[, hqc_idx]
meta <- meta[hqc_idx, ]

# keep the major tissue samples
keep_names <- table(meta$SMTS) %>% .[.>500] %>% names
mts_idx <- meta$SMTS %in% keep_names
meta <- meta[mts_idx, ]
gtpm <- gtpm[, mts_idx]

# keep the genes with high variance (top 75%)
var_rows <- mclapply(1:nrow(gtpm),
                     function(i) sd(gtpm[i,]),
                     mc.cores = 10) %>% unlist
hvg_idx <- (var_rows > quantile(var_rows, 0.75))
gtpm <- gtpm[hvg_idx, ]

## Obtain the embedding with UMAP
gexpr <- t(gtpm)
set.seed(123)
pca_res <- RSpectra::svds(gexpr, k = 30)
pcs_res <- pca_res$u %*% diag(pca_res$d)
emb <- uwot::umap(pcs_res,
                  n_components = 2,
                  n_neighbors=15,
                  min_dist = 0.5)
emb <- scale(emb, center = TRUE, scale = FALSE)
# df <- data.frame(
#   x = emb[,1],
#   y = emb[,2],
#   smts = meta$SMTS
# )
# ggplot(df, aes(x=x, y=y, col=smts)) + geom_point(size=1, show.legend = FALSE) + dark_theme_gray()

################################################################################
# Replace the tissue subcategory names
ts_names <- meta$SMTS %>% table %>% names
subts_names <- c()
for(tn in ts_names){
  add_names <-
    meta %>%
    filter(SMTS == tn) %>%
    pull(SMTSD) %>%
    table %>%
    names
  subts_names <- c(subts_names, sort(add_names))
}

new_subts_names <-
  c("Subcutaneous",
    "Visceral",
    "EBV-transformed lymphocytes",
    "Whole Blood",
    "Aorta",
    "Coronary",
    "Tibial",
    "Amygdala",
    "Anterior cingulate cortex",
    "Caudate",
    "Cerebellar Hemisphere",
    "Cerebellum",
    "Cortex",
    "Frontal Cortex",
    "Hippocampus",
    "Hypothalamus",
    "Nucleus accumbens",
    "Putamen",
    "Spinal cord",
    "Substantia nigra",
    "Gastroesophageal Junction",
    "Mucosa",
    "Muscularis",
    "Muscle - Skeletal",
    "Cultured fibroblasts",
    "Not Sun Exposed",
    "Sun Exposed")

subts_labels <- meta$SMTSD
for(i in seq_along(subts_names)){
  idx <- which(subts_labels == subts_names[i])
  subts_labels[idx] <- new_subts_names[i]
}
subts_names_levels <-
  c("Subcutaneous",
    "Visceral",
    "EBV-transformed lymphocytes",
    "Whole Blood",
    "Aorta",
    "Coronary",
    "Tibial",
    "Amygdala",
    "Anterior cingulate cortex",
    "Caudate",
    "Cerebellar Hemisphere",
    "Cerebellum",
    "Cortex",
    "Frontal Cortex",
    "Hippocampus",
    "Hypothalamus",
    "Nucleus accumbens",
    "Putamen",
    "Spinal cord",
    "Substantia nigra",
    "Gastroesophageal Junction",
    "Muscularis",
    "Mucosa",
    "Muscle - Skeletal",
    "Cultured fibroblasts",
    "Not Sun Exposed",
    "Sun Exposed")
meta$SMTSD <- factor(subts_labels, levels = subts_names_levels, ordered = FALSE)

## assign tissue colors and subcategory colors
tissue_names <- names(table(meta$SMTS))
tcolor <- ggColorHue(7)
names(tcolor) <- tissue_names
dmatrix <- table(meta$SMTS, meta$SMTSD)
tscolor <- c()
for(tn in tissue_names){
  ts <- names(which(dmatrix[tn,]!=0))
  ts <- ts[order(which(ts %in% new_subts_names))]
  if(length(ts) != 1){
    cur_tcolor <- tcolor[tn]
    color_base <- colorRampPalette(c("white", cur_tcolor, "black"))(100)[10:80]
    cut_color <- quantile(1:length(color_base),
                          probs = seq(0, 1, length.out = length(ts) + 2))
    cur_tscolor <- color_base[cut_color[c(-1,-length(cut_color))]]
  } else {
    cur_tscolor <- tcolor[tn]
  }
  names(cur_tscolor) <- ts
  tscolor <- c(tscolor, cur_tscolor)
}
tcol <- list(tcolor = tcolor,
             tscolor = tscolor)

## reorder the samples by the tissue subcategories
ind <- order(meta$SMTSD)
gtpm <- gtpm[, ind]
emb  <- emb[ind,]
meta <- meta[ind, ]

## save files
saveRDS(t(gtpm), file = "./GTEx/results/GTEx_expr.rds")
saveRDS(emb, file = "./GTEx/results/GTEx_emb.rds")
saveRDS(meta, file = "./GTEx/results/GTEx_meta.rds")
saveRDS(tcol, file = "./GTEx/results/GTEx_tcol.rds")

## Download GO terms (functional gene sets)
out = "./GTEx/data/GO_Biological_Process.txt"
url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023"
download.file(url, destfile = out)
