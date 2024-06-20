source("./R/functions.R")
library(fgsea)

## Input
emb <- readRDS(file="./GTEx/results/GTEx_emb.rds") # embedding matrix n*2
exp <- readRDS(file="./GTEx/results/GTEx_expr.rds") # expression matrix n*p
gsets <- fgsea::gmtPathways("./GTEx/data/GO_Biological_Process_2023.txt")

## Step1: filter and map the genes
inputs <- list(exp = exp, emb = emb)
mapped_res <- MapGenes(inputs)

## Step2: identify the enriched function
inputs <- mapped_res$inputs_enrich_test
set.seed(2024) # exclude the randomness introduced by the permutation test
enrich_res <- EnrichTest(inputs, gsets)

saveRDS(mapped_res, file = "./GTEx/results/mapped_genes.rds")
saveRDS(enrich_res, file = "./GTEx/results/enriched_terms.rds")
