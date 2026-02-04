library(scDesign3)
library(SingleCellExperiment)
library(ggplot2)
library(DuoClustering2018)
library(scran)
library(tidyverse)
theme_set(theme_bw())
############################
## Read in Reference Data ##
############################

sce <- get("sce_filteredExpr10_Zhengmix4eq")(metadata = FALSE)
colData(sce)$cell_type = as.factor(colData(sce)$phenoid)



ngene <- 200
logcounts(sce) <- log1p(counts(sce))
temp_sce <- modelGeneVar(sce)
chosen <- getTopHVGs(temp_sce, n = ngene)
sce <- sce[chosen,]


###############################################
## Simulation With Original Cell-Type Labels ##
###############################################


set.seed(123)
example_simu <- scdesign3(
  sce = sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  mu_formula = "cell_type",
  sigma_formula = "cell_type",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE,
  corr_formula = "cell_type",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = FALSE,
  parallelization = "pbmcmapply"
)

logcounts(sce) <- log1p(counts(sce))
simu_sce <- SingleCellExperiment(list(counts = example_simu$new_count), colData = example_simu$new_covariate)
logcounts(simu_sce) <- log1p(counts(simu_sce))

####################################################
## Simulation with modified cell-type proportions ##
####################################################

# 1

example_data <- construct_data(
  sce = sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  corr_by = "cell_type"
)

# 2

example_marginal <- fit_marginal(
  data = example_data,
  predictor = "gene",
  mu_formula = "cell_type",
  sigma_formula = "cell_type",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE,
  parallelization = "pbmcmapply"
)

# 3 

set.seed(123)
example_copula <- fit_copula(
  sce = sce,
  assay_use = "counts",
  marginal_list = example_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 2,
  input_data = example_data$dat
)

# 4

table(colData(sce)$cell_type) / dim(sce)[2]

ct_prop <- c(0, 0, 0.2,0.8)
unique_combined <- example_data$dat %>% tidyr::expand(nesting(cell_type, corr_group))
new_ct <- as.data.frame(lapply(unique_combined, rep,round(ct_prop*dim(sce)[2])))
head(new_ct)
table(new_ct$cell_type)/dim(new_ct)[1]

example_para <- extract_para(
  sce = sce,
  marginal_list = example_marginal,
  n_cores = 1,
  family_use = "nb",
  new_covariate = new_ct,
  data = example_data$dat
)

# 5

set.seed(123)
example_newcount <- simu_new(
  sce = sce,
  mean_mat = example_para$mean_mat,
  sigma_mat = example_para$sigma_mat,
  zero_mat = example_para$zero_mat,
  quantile_mat = NULL,
  copula_list = example_copula$copula_list,
  n_cores = 1,
  family_use = "nb",
  input_data = example_data$dat,
  new_covariate = new_ct,
  important_feature = example_copula$important_feature,
  filtered_gene = example_data$filtered_gene
)

logcounts(sce) <- log1p(counts(sce))
simu_sce2 <- SingleCellExperiment(list(counts = example_newcount), colData = data.frame(cell_type = new_ct$cell_type))
logcounts(simu_sce2) <- log1p(counts(simu_sce2))

###################
## Visualization ##
###################

set.seed(123)
compare_figure <- plot_reduceddim(ref_sce = sce, 
                                  sce_list = list(simu_sce, simu_sce2), 
                                  name_vec = c("Reference", "Same cell-type proportions", "Modified cell-type proportions"),
                                  assay_use = "logcounts", 
                                  if_plot = TRUE, 
                                  color_by = "cell_type", 
                                  n_pc = 20)
plot(compare_figure$p_umap)

sessionInfo()
