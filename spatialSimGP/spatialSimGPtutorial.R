library(MASS)
library(SpatialExperiment)
#BiocManager::install("STexampleData")
library(STexampleData)
library(ggplot2)
library(spatialSimGP)

spe_demo <- Visium_mouseCoronal()

colData(spe_demo)$subset <- ifelse(colData(spe_demo)$array_row > 20 & colData(spe_demo)$array_row < 65 & colData(spe_demo)$array_col > 30 & colData(spe_demo)$array_col < 65, TRUE, FALSE)
spe_demo <- spe_demo[, colData(spe_demo)$subset]
coords <- spatialCoords(spe_demo)

n_genes <- 5
proportion <- 0.4
range_sigma.sq <- c(1.5, 3)
range_beta <- c(3, 7)

length_scale <- 60

set.seed(16)
spe <- spatial_simulate(n_genes, proportion, coords, range_sigma.sq, range_beta,length_scale, length_scale_option = "fixed")

df <- as.data.frame(cbind(spatialCoords(spe), expr = counts(spe)[1, ]))

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = expr)) + geom_point(size = 2.2) + coord_fixed() + scale_y_reverse() + scale_color_gradient(low = "gray90", high = "blue", trans = "sqrt", breaks = range(df$expr), name = "counts") + theme_bw() + theme(plot.title = element_text(face = "italic"), panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

length_scale <- c(60, 40, 20, 80, 100)

set.seed(1)
spe <- spatial_simulate(n_genes, proportion, coords,
                        range_sigma.sq, range_beta,
                        length_scale, length_scale_option = "unique")


n_side <- 20

# x and y coordinates for the grid
x_coords <- rep(1:n_side, each = n_side)
y_coords <- rep(1:n_side, times = n_side)

# combine into a matrix
coords <- cbind(x_coords, y_coords)
colnames(coords) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")

# run the simulation
set.seed(1)
length_scale <- 60
spe <- spatial_simulate(n_genes, proportion, coords,
                        range_sigma.sq, range_beta,
                        length_scale, length_scale_option = "fixed")


df <- as.data.frame(cbind(spatialCoords(spe), expr = counts(spe)[1, ]))

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = expr)) + 
  geom_point(size = 5) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "blue", 
                       trans = "sqrt", breaks = range(df$expr), 
                       name = "counts") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "italic"), 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())


set.seed(123) 
n_genes <- 1 
proportion <- 0 
range_sigma.sq <- c(1, 1)
range_beta <- c(3, 3)
length_scale <- 60

spe <- spatial_simulate(n_genes, proportion, coords,
                        range_sigma.sq, range_beta,
                        length_scale, length_scale_option = "fixed")

sessionInfo()
