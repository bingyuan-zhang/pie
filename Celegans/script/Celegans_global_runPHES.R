source("./R/functions.R")

## Prepare inputs
emb <- readRDS("./Celegans/results/global_emb.rds")
exp <- readRDS("./Celegans/results/global_exp.rds")
gsets <- fgsea::gmtPathways("./Celegans/data/GO_Biological_Process_2018.txt")

## Step1: filter and map the genes
inputs <- list(exp = exp, emb = emb)
mapped <- MapGenes(inputs)

## Step2: identify the enriched function
inputs <- mapped$inputs_enrich_test
set.seed(2024)
enriched <- EnrichTest(inputs, gsets)

## Save
saveRDS(enriched, "./Celegans/results/global_enriched_terms.rds")
saveRDS(mapped, "./Celegans/results/global_mapped_genes.rds")

