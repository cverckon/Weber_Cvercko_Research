
#########################################
### varying sdmax in GenerateTissue() ###
#########################################

library(scSpatialSIM)
set.seed(333) #reproducibility

custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

# create objs
sim_obj_1 <- CreateSimulationObject(sims = 9, cell_types = 1, window = custom_window)
sim_obj_2 <- CreateSimulationObject(sims = 9, cell_types = 1, window = custom_window)
sim_obj_3 <- CreateSimulationObject(sims = 9, cell_types = 1, window = custom_window)
sim_obj_4 <- CreateSimulationObject(sims = 9, cell_types = 1, window = custom_window)
sim_obj_5 <- CreateSimulationObject(sims = 9, cell_types = 1, window = custom_window)

summary(sim_obj_1)
summary(sim_obj_2)
summary(sim_obj_3)
summary(sim_obj_4)
summary(sim_obj_5)


# generate spatial patterns
sim_obj_1 <- GenerateSpatialPattern(sim_obj_1, 25)
sim_obj_2 <- GenerateSpatialPattern(sim_obj_2, 25)
sim_obj_3 <- GenerateSpatialPattern(sim_obj_3, 25)
sim_obj_4 <- GenerateSpatialPattern(sim_obj_4, 25)
sim_obj_5 <- GenerateSpatialPattern(sim_obj_5, 25)

summary(sim_obj_1)
summary(sim_obj_2)
summary(sim_obj_3)
summary(sim_obj_4)
summary(sim_obj_5)

plot(sim_obj_1, what = "Patterns", ncol = 1, nrow = 1, which = 1) 
plot(sim_obj_2, what = "Patterns", ncol = 1, nrow = 1, which = 1) 
plot(sim_obj_3, what = "Patterns", ncol = 1, nrow = 1, which = 1) 
plot(sim_obj_4, what = "Patterns", ncol = 1, nrow = 1, which = 1)
plot(sim_obj_5, what = "Patterns", ncol = 1, nrow = 1, which = 1)


# generate tissues
sim_obj_1 <- GenerateTissue(sim_obj_1, density_heatmap = T, sdmax = 2, step_size = 0.1, cores = 1)
sim_obj_2 <- GenerateTissue(sim_obj_2, density_heatmap = T, sdmax = 1, step_size = 0.1, cores = 1)
sim_obj_3 <- GenerateTissue(sim_obj_3, density_heatmap = T, sdmax = 0.5, step_size = 0.1, cores = 1)
sim_obj_4 <- GenerateTissue(sim_obj_4, density_heatmap = T, sdmax = 0.25, step_size = 0.1, cores = 1)
sim_obj_5 <- GenerateTissue(sim_obj_5, density_heatmap = T, sdmax = 0.1, step_size = 0.1, cores = 1)


summary(sim_obj_1)
summary(sim_obj_2)
summary(sim_obj_3)
summary(sim_obj_4)
summary(sim_obj_5)


PlotSimulation(sim_obj_1, which = 1:4, ncol = 2, nrow = 2, what = "tissue heatmap")
PlotSimulation(sim_obj_2, which = 1:4, ncol = 2, nrow = 2, what = "tissue heatmap")
PlotSimulation(sim_obj_3, which = 1:4, ncol = 2, nrow = 2, what = "tissue heatmap")
PlotSimulation(sim_obj_4, which = 1:4, ncol = 2, nrow = 2, what = "tissue heatmap")
PlotSimulation(sim_obj_5, which = 1:4, ncol = 2, nrow = 2, what = "tissue heatmap")
warnings()


# final without generating holes
sim_obj_1 <- GenerateCellPositivity(sim_obj_1, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)
sim_obj_2 <- GenerateCellPositivity(sim_obj_2, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)
sim_obj_3 <- GenerateCellPositivity(sim_obj_3, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)
sim_obj_4 <- GenerateCellPositivity(sim_obj_4, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)


PlotSimulation(sim_obj_1, which = 1, what = "whole core")
PlotSimulation(sim_obj_2, which = 1, what = "whole core")
PlotSimulation(sim_obj_3, which = 1, what = "whole core")
PlotSimulation(sim_obj_4, which = 1, what = "whole core")


# generating holes
sim_obj_1 <- GenerateHoles(sim_obj_1, density_heatmap = T, step_size = 0.1, cores = 1)
sim_obj_2 <- GenerateHoles(sim_obj_2, density_heatmap = T, step_size = 0.1, cores = 1)
sim_obj_3 <- GenerateHoles(sim_obj_3, density_heatmap = T, step_size = 0.1, cores = 1)
sim_obj_4 <- GenerateHoles(sim_obj_4, density_heatmap = T, step_size = 0.1, cores = 1)


summary(sim_obj_1)
summary(sim_obj_2)
summary(sim_obj_3)
summary(sim_obj_4)

PlotSimulation(sim_obj_1, which = 1:8, ncol = 2, nrow = 2, what = "hole heatmap")
PlotSimulation(sim_obj_2, which = 1:8, ncol = 2, nrow = 2, what = "hole heatmap")
PlotSimulation(sim_obj_3, which = 1:8, ncol = 2, nrow = 2, what = "hole heatmap")
PlotSimulation(sim_obj_4, which = 1:8, ncol = 2, nrow = 2, what = "hole heatmap")


# Generate finals
sim_obj_1 <- GenerateCellPositivity(sim_obj_1, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)
sim_obj_1 <- GenerateCellPositivity(sim_obj_2, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)
sim_obj_1 <- GenerateCellPositivity(sim_obj_3, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)
sim_obj_1 <- GenerateCellPositivity(sim_obj_4, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)


summary(sim_obj_1)
summary(sim_obj_2)
summary(sim_obj_3)
summary(sim_obj_4)

PlotSimulation(sim_obj_1, which = 1, what = "whole core")
PlotSimulation(sim_obj_2, which = 1, what = "whole core")
PlotSimulation(sim_obj_3, which = 1, what = "whole core")
PlotSimulation(sim_obj_4, which = 1, what = "whole core")