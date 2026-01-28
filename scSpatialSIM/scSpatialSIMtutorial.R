if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("FridleyLab/mIFsim@v0.1.3.4")

library(scSpatialSIM)
set.seed(333) #reproducibility

custom_window = spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

sim_object = CreateSimulationObject(sims = 9, cell_types = 1, window = custom_window)

summary(sim_object)


sim_object = GenerateSpatialPattern(sim_object)
summary(sim_object)
plot(sim_object, what = "Patterns", ncol = 1, nrow = 1, which = 1)#print only first point pattern

sim_object = GenerateTissue(sim_object, density_heatmap = T, step_size = 0.1, cores = 1)
summary(sim_object)


PlotSimulation(sim_object, which = 1:4, ncol = 2, nrow = 2, what = "tissue heatmap")

sim_object = GenerateHoles(sim_object, density_heatmap = T, step_size = 0.1, cores = 1)
summary(sim_object)

PlotSimulation(sim_object, which = 1:8, ncol = 2, nrow = 2, what = "hole heatmap")


sim_object = GenerateCellPositivity(sim_object, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)
summary(sim_object)

PlotSimulation(sim_object, which = 1, what = "whole core")


##########################
## Bivariate Simulation ##
##########################

#set seed
set.seed(333)
#create the new object
bivariate_sim = CreateSimulationObject(sims = 5, cell_types = 2) %>%
  #produce the point pattern
  GenerateSpatialPattern() %>%
  #make tissues
  GenerateTissue(density_heatmap = T, step_size = 0.1, cores = 1)

bivariate_sim_tmp  = GenerateCellPositivity(bivariate_sim, k = 4,
                                            sdmin = 3, sdmax = 5,
                                            density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                            shift = 0)

PlotSimulation(bivariate_sim_tmp, which = 1, what = "whole core")

bivariate_sim_tmp  = GenerateCellPositivity(bivariate_sim, k = 4,
                                            sdmin = 3, sdmax = 5,
                                            density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                            shift = 1)

PlotSimulation(bivariate_sim_tmp, which = 1, what = "whole core")


####################
## Exporting Data ##
####################

spatial_list = CreateSpatialList(sim_object = bivariate_sim_tmp)
head(spatial_list[[1]])

single_dataframe = CreateSpatialList(sim_object = bivariate_sim_tmp, single_df = TRUE)
head(single_dataframe)

summary_data = SummariseSpatial(spatial_list = spatial_list, markers = c("Cell 1 Assignment", "Cell 2 Assignment"))
head(summary_data)

