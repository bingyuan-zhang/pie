library(dplyr)

source("./R/functions.R")

emb <- readRDS("./Celegans/results/global_emb.rds")
exp <- readRDS("./Celegans/results/global_exp.rds")
meta <- readRDS("./Celegans/results/global_meta.rds")

## focus on the subset
cell_type <- meta$cell.type
focus <- c("ASJ", "AUA", "ASE_parent", "ASE", "ASEL", "Neuroblast_ASJ_AUA", "Neuroblast_ASE_ASJ_AUA")
ind1 <- cell_type %in% focus # condition 1
ind2 <- emb[,1] < -3.5 # condition 2

## embedding structure
emb <- emb[ind1 & ind2, ] %>% scale(center = TRUE, scale = FALSE)
exp <- exp[ind1 & ind2, ]
meta <- meta[ind1 & ind2, ]

## computational pipeline
inputs <- list(exp = exp, emb = emb)
gsets <- fgsea::gmtPathways("./Celegans/data/GO_Biological_Process_2018.txt")
mapped <- MapGenes(inputs)
set.seed(2024)
enriched <- EnrichTest(inputs = mapped$inputs_enrich_test, gene_sets = gsets)

## save data
saveRDS(emb, file = "./Celegans/results/local_emb.rds")
saveRDS(exp, file = "./Celegans/results/local_exp.rds")
saveRDS(meta, file = "./Celegans/results/local_meta.rds")
saveRDS(enriched, file = "./Celegans/results/local_enriched_terms.rds")
saveRDS(mapped, file = "./Celegans/results/local_mapped_genes.rds")







