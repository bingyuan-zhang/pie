library(readr)
library(dplyr)
library(stringr)

source("./Rscript/functions.R")
exp    <- readRDS("./GTEx/results/GTEx_expr.rds")
terms  <- readRDS("./GTEx/results/enriched_terms.rds")
mapped <- readRDS("./GTEx/results/mapped_genes.rds")
term_by_intv <- readRDS("./GTEx/results//term_intv.rds")
gsets <- fgsea::gmtPathways("./GTEx/data/GO_Biological_Process_2023.txt")

## Enriched GO terms
terms_pd    <- terms$pd # projected direction (coordinates) of the eigengene
terms_pval  <- terms$pval # p-value of the GSEA test
terms_padj  <- terms$padj # adjusted p-value of the GSEA test
terms_egene <- terms$egene # eigengene vector with length n
terms_name  <- terms$name %>% # name of the GO term
  strsplit(split = "\\(") %>%
  map_vec(~ head(.x, 1)) %>%
  str_remove("\\)") %>%
  str_remove(" $")

## tissue specific genes
tissue_category <- read_csv("./GTEx/data/tissue_category.csv")
is_ts_enriched <- tissue_category$`RNA tissue specificity` == "Tissue enriched"
is_gr_enriched <- tissue_category$`RNA tissue specificity` == "Group enriched"
is_ts_enhanced <- tissue_category$`RNA tissue specificity` == "Tissue enhanced"
sub <- tissue_category %>% filter(is_ts_enhanced |  is_ts_enriched | is_gr_enriched)

sub_skin  <-
  "skin" %>% grep(., sub$`RNA tissue specific nTPM`) %>% sub[., ]
sub_brain <-
  "brain" %>% grep(., sub$`RNA tissue specific nTPM`) %>% sub[., ]
sub_esophagus <-
  "esophagus" %>% grep(., sub$`RNA tissue specific nTPM`) %>% sub[., ]
sub_skeletal <-
  "skeletal" %>% grep(., sub$`RNA tissue specific nTPM`) %>% sub[., ]
sub_adipose <-
  "adipose" %>% grep(., sub$`RNA tissue specific nTPM`) %>% sub[., ]

glist <- list(
  "Skin" = intersect(sub_skin$Gene, colnames(exp)),
  "Brain" = intersect(sub_brain$Gene, colnames(exp)),
  "Esophagus" = intersect(sub_esophagus$Gene, colnames(exp)),
  "Muscle" = intersect(sub_skeletal$Gene, colnames(exp)),
  "Adipose" = intersect(sub_adipose$Gene, colnames(exp))
)

# prepare the selected genes
true_g = mapped$select
# prepare the enriched functions
true_f  = terms$name

## figure 5a
ov <- lapply(glist, function(x) sum(x%in% mapped$select)) %>% unlist
## figure 5a
ov <- lapply(glist, function(x) sum(x %in% mapped$select)) %>% unlist  # figure 5b still uses ov
# a gene can be specific to several tissues, so count unique genes
n_ts  <- sum(unique(unlist(glist)) %in% mapped$select)
n_nts <- length(mapped$select) - n_ts
n_ts; n_nts
data <- data.frame(
  tsgenes = c("Tissue-specific", "Non-tissue-specific"),
  number  = c(n_ts, n_nts),
  color   = c("#F8766D", "grey90")
)

data <- data %>%
  mutate(tsgenes = factor(tsgenes, levels = tsgenes))

figure_5a <-
  ggplot(data, aes(x = "", y = number, fill = tsgenes)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = number),
            position = position_stack(vjust = 0.5),
            size = 18/2.845,
            col = "black") +
  scale_fill_manual(values = data$color) +
  labs(fill = "PIE's mapped genes ", x = "", y = "") +
  theme_void() +
  theme(text = element_text(family = "Arial"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))

figure_5a

## figure_5b
find.pvals <- function(hpa){
  pie <- mapped$select
  background <- colnames(exp)
  hpa_cond <- background %in% hpa
  pie_cond <- background %in% pie

  a <- (hpa_cond & pie_cond) %>% sum
  b <- (hpa_cond & !pie_cond) %>% sum()
  c <- (!hpa_cond & pie_cond) %>% sum()
  d <- (!hpa_cond & !pie_cond) %>% sum()

  test_mat <- matrix(c(a,b,c,d), nrow=2)
  colnames(test_mat) <- c("is.hpa", "isnot.hpa")
  rownames(test_mat) <- c("is.pie.mapped", "isnot.pie.mapped")
  fisher_res <- fisher.test(test_mat)
  fisher_res$p.value
}

pval <- lapply(glist, function(x) find.pvals(x)) %>% unlist
df <- data.frame(ts = names(pval), pval = pval, ov = ov)
df$log_pval <- -log10(df$pval)
df <- df %>%
  mutate(ts = reorder(ts, log_pval, FUN = function(x)-x),
         pval_text = sprintf("%.2g", pval),
         col = c("#FB61D7", "#00C094", "#00B6EB", "#A58AFF", "#F8766D"))
figure_5b <-
  ggplot(df, aes(x = ts,
                 y = log_pval,
                 fill = col)) +
  geom_bar(stat = "identity",
           width = 0.7) +
  geom_text(aes(label = pval_text),
            vjust = -0.2,
            color = "black",
            size = 20/2.845) +
  geom_text(x = 5,
            y = -log10(0.05),
            label = "p = 0.05",
            vjust = -1,
            size = 20/2.845,
            fontface = "bold",
            color = "red") +
  geom_hline(yintercept = -log10(0.05),
             linetype = 6,
             color = "red",
             size = 1) +
  scale_fill_identity() +
  labs(x = "", y = "-log10(p-value)", fill = "") +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

figure_5b


## figure 5c
library(enrichR)
setEnrichrSite("Enrichr")
enriched = map(glist, ~ enrichr(.x, "GO_Biological_Process_2023"))
range = 1:20 # compare the top 20
ovlp <- matrix(0, nrow=length(range), ncol=5)
colnames(ovlp) <- names(glist)
enriched$Skin$GO_Biological_Process_2023$Term[1:20]
for(i in 1:5){
  A <- enriched[[i]]$GO_Biological_Process_2023$Term[range] %>%
    strsplit(split = "\\(") %>%
    map_vec(~ head(.x, 1)) %>%
    str_remove("\\)") %>%
    str_remove(" $")
  ovlp[,i] <- ifelse(A %in% terms_name, "Yes", "No")
}

df_ovlp_func <- as.data.frame.table(ovlp)
colnames(df_ovlp_func) <- c("top", "tissue", "ovlp")
df_ovlp_func$top <- as.numeric(df_ovlp_func$top)
df_ovlp_func$tissue <-
  factor(df_ovlp_func$tissue,
         levels = c("Adipose", "Skin", "Muscle", "Brain", "Esophagus"),
         ordered = FALSE)
figure_5c <-
  ggplot(df_ovlp_func, aes(x = top, y = tissue, fill = ovlp)) +
  geom_tile(color = "white", width = 0.8, height = 0.8) +
  scale_fill_manual(values = c("Yes" = "red2", "No" = "grey80")) +
  coord_fixed(ratio = 1) +
  labs(y = "",
       x = "Overlap of the top tissue specific functions",
       fill = "PIE's mapped functions") +
  geom_text(data = df_ovlp_func[df_ovlp_func$top==1, ],
            aes(y = tissue,
                x = 0,
                label = tissue),
            size = 18/2.845,
            hjust = 1,
            vjust = 0.5,
            col = "black") +
  geom_text(data = df_ovlp_func[df_ovlp_func$tissue=="Skin",],
            aes(y = 0,
                x = top,
                label = top),
            size = 15/2.845,
            hjust = 0.5,
            vjust = -0.3,
            col = "black") +
  xlim(c(-3,21)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(text = element_text(family = "Arial"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 18, vjust = -0.5),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

figure_5c


## Evaluate Reproducibility
niter <- 10
spsample_vec <- c(0.2, 0.4, 0.6, 0.8)
result_gmat <- array(dim=c(4, 10))
result_fmat <- array(dim=c(4, 10))
gexpr <- readRDS(, file = "./GTEx/results/GTEx_expr.rds")

split_sample <- function(gexpr, prop) {
  ind <- sample(nrow(gexpr), size = prop* nrow(gexpr), replace = FALSE)
  # submatrix
  gexpr_new <- gexpr[ind, ]
  pca_res <- RSpectra::svds(gexpr_new, k = 30)
  pcs_res <- pca_res$u %*% diag(pca_res$d)
  
  # obtain new UMAP embedding
  res_umap <- uwot::umap(pcs_res,
                         n_components = 2,
                         n_neighbors = 15,
                         min_dist = 0.5,
                         seed = 123)
  emb_umap_new <- scale(res_umap, center = TRUE, scale = FALSE)
  list(exp = gexpr_new, emb = emb_umap_new)
}

set.seed(2024)
for(j in seq_along(spsample_vec)) {
  for (iter in 1:niter) {
    print(paste0(j, " ", iter))
    prop = spsample_vec[j]
    split_inputs = split_sample(gexpr, prop)

    mapped = MapGenes(split_inputs)

    enriched = EnrichTest(mapped$inputs_enrich_test, gsets)

    select_g = mapped$select
    select_f = enriched$name
    recover_g_rate = length(intersect(select_g, true_g)) / length(true_g)
    recover_f_rate = length(intersect(select_f, true_f)) / length(true_f)

    result_gmat[j, iter] = recover_g_rate
    result_fmat[j, iter] = recover_f_rate
  }
}
saveRDS(result_gmat, file="./GTEx/results/eval_g.rds")
saveRDS(result_fmat, file="./GTEx/results/eval_f.rds")

## figure 5d
result_gmat <- readRDS("./GTEx/results/eval_g.rds")
result_fmat <- readRDS("./GTEx/results/eval_f.rds")

indices <- which(result_gmat == result_gmat, arr.ind = TRUE)
df_g <- data.frame(sample_prop = c("20%", "40%", "60%", "80%")[indices[, 1]],
                   value = result_gmat[result_gmat == result_gmat])
figure_5d <-
  ggplot(df_g, aes(y = value, x = sample_prop)) +
  geom_boxplot(show.legend = FALSE,
               outlier.shape = NA,
               width = 0.6,
               fill = "yellow") +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.6,
    color = "red",
    size = 1
  ) +
  geom_point(
    position = position_jitter(width = 0.2),
    size = 3,
    color = "black",
    alpha = 1
  ) +
  labs(y = "Reproducibility ratio of genes",
       x = "Sample proportion") +
  ylim(c(0.8, 1)) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20,
                                color = "black",
                                margin = margin(r = 15)),
    axis.title.x = element_text(size = 20,
                                color = "black",
                                margin = margin(t = 15)),

    axis.ticks = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    ),
    panel.background = element_blank()
  )
figure_5d

## figure 5e
indices <- which(result_fmat == result_fmat, arr.ind = TRUE)
df_f <- data.frame(
  sample_prop = c("20%", "40%", "60%", "80%")[indices[, 1]],
  value = result_fmat[result_fmat == result_fmat]
)
figure_5e <-
  ggplot(df_f, aes(y=value, x=sample_prop)) +
  geom_boxplot(show.legend = FALSE,
               outlier.shape = NA,
               width = 0.6,
               fill = "yellow") +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.6,
    color = "red",
    size = 1
  ) +
  geom_point(
    position = position_jitter(width = 0.2),
    size = 3,
    color = "black",
    alpha = 1
  ) +
  labs(y = "Reproducibility ratio of functions",
       x = "Sample proportion") +
  ylim(c(0.5, 1)) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20,
                                color = "black",
                                margin = margin(r = 15)),
    axis.title.x = element_text(size = 20,
                                color = "black",
                                margin = margin(t = 15)),

    axis.ticks = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    ),
    panel.background = element_blank()
  )

figure_5e


ggsave("./GTEx/figures/Figure_5A.svg",
       plot = figure_5a, width = 5.0, height = 3.5,
       units = "in", device = svglite::svglite, bg = "white")

ggsave("./GTEx/figures/Figure_5B.svg",
       plot = figure_5b, width = 4.5, height = 3.5,
       units = "in", device = svglite::svglite, bg = "white")

ggsave("./GTEx/figures/Figure_5C.svg",
       plot = figure_5c, width = 9.0, height = 2.8,
       units = "in", device = svglite::svglite, bg = "white")

ggsave("./GTEx/figures/Figure_5D.svg",
       plot = figure_5d, width = 4.5, height = 4.0,
       units = "in", device = svglite::svglite, bg = "white")

ggsave("./GTEx/figures/Figure_5E.svg",
       plot = figure_5e, width = 4.5, height = 4.0,
       units = "in", device = svglite::svglite, bg = "white")
