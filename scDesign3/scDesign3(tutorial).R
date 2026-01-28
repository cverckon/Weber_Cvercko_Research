library(scDesign3)
library(ggplot2)
library(BiocManager)
library(SingleCellExperiment)
library(dplyr)
library(viridis)


theme_set(theme_bw())

example_sce <- readRDS((url("https://figshare.com/ndownloader/files/40582019")))


readRDS(url("https://figshare.com/ndownloader/files/40582019"))

example_sce <- readRDS('VISIUM_sce.rds')
print(example_sce)

example_sce <- example_sce[1:10, ];example_sce


set.seed(123)
example_simu <- scdesign3(
  sce = example_sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  mu_formula = "s(spatial1, spatial2, bs = 'gp', k= 400)",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = FALSE,
  parallelization = "pbmcapply"
)

simu_sce <- SingleCellExperiment(list(counts = example_simu$new_count), colData = example_simu$new_covariate)
logcounts(simu_sce) <- log1p(counts(simu_sce))


VISIUM_dat_test <- data.frame(t(log1p(counts(example_sce)))) %>% as_tibble() %>% dplyr::mutate(X = colData(example_sce)$spatial1, Y = colData(example_sce)$spatial2) %>% tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "Reference")
VISIUM_dat_scDesign3 <- data.frame(t(log1p(counts(simu_sce)))) %>% as_tibble() %>% dplyr::mutate(X = colData(simu_sce)$spatial1, Y = colData(simu_sce)$spatial2) %>% tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "scDesign3")
VISIUM_dat <- bind_rows(VISIUM_dat_test, VISIUM_dat_scDesign3) %>% dplyr::mutate(Method = factor(Method, levels = c("Reference", "scDesign3")))

VISIUM_dat %>% filter(Gene %in% rownames(example_sce)[1:5]) %>% ggplot(aes(x = X, y = Y, color = Expression)) + geom_point(size = 0.5) + scale_colour_gradientn(colors = viridis_pal(option = "magma")(10), limits=c(0, 8)) + coord_fixed(ratio = 1) + facet_grid(Method ~ Gene )+ theme_gray()


