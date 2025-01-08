source("./Rscript/functions.R")

## Input
emb <- readRDS(file="./GTEx/results/GTEx_emb.rds") # embedding matrix n*2
exp <- readRDS(file="./GTEx/results/GTEx_expr.rds") # expression matrix n*p
gene_sets <- fgsea::gmtPathways("./GTEx/data/GO_Biological_Process_2023.txt")

# Compare the runtime in PIE with fixed n and p.
runtime <- c("Step1.scmarker"=0, "Step1.map"=0,
             "Step2.svd"=0, "Step2.map"=0, "Step2.gsea"=0)

for(rep in 1:3) {
  ## Step1: filter and map the genes
  inputs <- list(exp = exp, emb = emb)
  mapped_res <- MapGenes(inputs)

  runtime["Step1.scmarker"] = runtime["Step1.scmarker"] + mapped_res$time_scmarker
  runtime["Step1.map"] = runtime["Step1.map"] + mapped_res$time_map

  ## Step2: identify the enriched function
  inputs <- mapped_res$inputs_enrich_test
  set.seed(2024) # exclude the randomness introduced by the permutation test
  enrich_res <- EnrichTest(inputs, gene_sets)

  runtime["Step2.svd"] = runtime["Step2.svd"] + enrich_res$time_svd
  runtime["Step2.map"] = runtime["Step2.map"] + enrich_res$time_map
  runtime["Step2.gsea"] = runtime["Step2.gsea"] + enrich_res$time_gsea
}
runtime_fixed_n_and_p <- runtime / 3
saveRDS(runtime_fixed_n_and_p, file="./GTEx/results/runtime_fixed_n_and_p.rds")

# Compare the runtime with changing n and fixed p
emb <- readRDS(file="./GTEx/results/GTEx_emb.rds") # embedding matrix n*2
exp <- readRDS(file="./GTEx/results/GTEx_expr.rds") # expression matrix n*p

runtime <- matrix(0, nrow=2, ncol=4)
colnames(runtime) <- c("1000", "2000", "3000", "4000")
rownames(runtime) <- c("step1", "step2")
sample_size <- c(1000, 2000, 3000, 4000)
for(i in 1:4) {
  sampled_ind <- sample(nrow(exp), sample_size[i])
  emb_cur <- emb[sampled_ind,]
  exp_cur <- exp[sampled_ind,]

  for(rep in 1:3) {
    ## Step1: filter and map the genes
    inputs <- list(exp = exp_cur, emb = emb_cur)
    tic <- proc.time()
    mapped_res <- MapGenes(inputs)
    toc <- proc.time()
    runtime[1,i] <- runtime[1,i] + (toc-tic)[3]

    ## Step2: identify the enriched function
    inputs <- mapped_res$inputs_enrich_test
    set.seed(2024) # exclude the randomness introduced by the permutation test
    tic <- proc.time()
    enrich_res <- EnrichTest(inputs, gene_sets)
    toc <- proc.time()
    runtime[2,i] <- runtime[2,i] + (toc-tic)[3]

  }
}
runtime_diff_n_fixed_p <- runtime / 3
saveRDS(runtime_diff_n_fixed_p, file="./GTEx/results/runtime_diff_n_fixed_p.rds")

# Compare the runtime with fixed n and different p
emb <- readRDS(file="./GTEx/results/GTEx_emb.rds") # embedding matrix n*2
exp <- readRDS(file="./GTEx/results/GTEx_expr.rds") # expression matrix n*p

runtime <- matrix(0, nrow=2, ncol=4)
colnames(runtime) <- c("200", "400", "600", "800")
rownames(runtime) <- c("step1", "step2")
feature_size <- c(200, 400, 600, 800)
for(i in 1:4) {
  sampled_ind <- sample(ncol(exp), feature_size[i])
  exp_cur <- exp[, sampled_ind]

  for(rep in 1:3) {
    ## Step1: filter and map the genes
    inputs <- list(exp = exp_cur, emb = emb)
    tic <- proc.time()
    mapped_res <- MapGenes(inputs)
    toc <- proc.time()
    runtime[1,i] <- runtime[1,i] + (toc-tic)[3]

    ## Step2: identify the enriched function
    inputs <- mapped_res$inputs_enrich_test
    set.seed(2024) # exclude the randomness introduced by the permutation test
    tic <- proc.time()
    enrich_res <- EnrichTest(inputs, gene_sets)
    toc <- proc.time()
    runtime[2,i] <- runtime[2,i] + (toc-tic)[3]

  }
}
runtime_fixed_n_diff_p <- runtime / 3
saveRDS(runtime_fixed_n_diff_p, file="./GTEx/results/runtime_fixed_n_diff_p.rds")


# plot the figures
library(ggplot2)
runtime_fixed_n_and_p <- readRDS(file="./GTEx/results/runtime_fixed_n_and_p.rds")
df <- data.frame(
  Step = factor(c( "Step1.scmarker", "Step1.map", "Step2.svd", "Step2.map", "Step2.gsea"),
                levels = c("Step1.scmarker", "Step1.map", "Step2.svd", "Step2.map", "Step2.gsea"))
)
df$Time <- runtime_fixed_n_and_p[as.character(df$Step)]
my_colors <- c(
  "Step1.map"      = "#000080",
  "Step1.scmarker" = "#1E90FF",
  "Step2.svd"      = "#FFB3B3",  # 浅红
  "Step2.map"      = "#FF6666",  # 中红
  "Step2.gsea"     = "#CC0000"   # 深红
)
ggplot(df, aes(x = Step, y = Time, fill = Step)) +
  geom_col() +
  scale_fill_manual(values = my_colors) +
  labs(x = NULL, y = "Second", title = "Runtime of steps") +
  theme_bw(base_size = 14)

runtime_fixed_n_diff_p <- readRDS(file="./GTEx/results/runtime_fixed_n_diff_p.rds")
df <- data.frame(
  feature_size = c(1000, 2000, 3000, 4000),
  Step1       = runtime_fixed_n_diff_p[1,],
  Step2       = runtime_fixed_n_diff_p[2,]
)
df_long <- tidyr::pivot_longer(
  data      = df,
  cols      = c("Step1", "Step2"),
  names_to  = "step",
  values_to = "value"
)
ggplot(df_long, aes(x = feature_size, y = value, color = step)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Step1" = "blue", "Step2" = "red")) +
  labs(
    x = "Number of Features",
    y = "Second",
    color = "Step",
    title = "Runtime with fixed samples"
  ) +
  theme_bw(base_size = 14)

runtime_diff_n_fixed_p <- readRDS(file="./GTEx/results/runtime_diff_n_fixed_p.rds")
df <- data.frame(
  sample_size = c(1000, 2000, 3000, 4000),
  Step1       = runtime_diff_n_fixed_p[1,],
  Step2       = runtime_diff_n_fixed_p[2,]
)
df_long <- tidyr::pivot_longer(
  data      = df,
  cols      = c("Step1", "Step2"),
  names_to  = "step",
  values_to = "value"
)
ggplot(df_long, aes(x = sample_size, y = value, color = step)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Step1" = "blue", "Step2" = "red")) +
  labs(
    x = "Number of Samples",
    y = "Second",
    color = "Step",
    title = "Runtime with fixed features"
  ) +
  theme_bw(base_size = 14)

