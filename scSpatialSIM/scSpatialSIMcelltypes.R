library(scSpatialSIM)
set.seed(333) #reproducibility

custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

# create obj
sim_obj_1 <- CreateSimulationObject(sims = 1, cell_types = 5, window = custom_window)
summary(sim_obj_1)


# generate spatial pattern
sim_obj_1 <- GenerateSpatialPattern(sim_obj_1, 25)
summary(sim_obj_1)

plot(sim_obj_1, what = "Patterns", ncol = 1, nrow = 1, which = 1) 


# generate tissue
sim_obj_1 <- GenerateTissue(sim_obj_1, density_heatmap = T, step_size = 0.1, cores = 1)
summary(sim_obj_1)


PlotSimulation(sim_obj_1, which = 1, ncol = 2, nrow = 2, what = "tissue heatmap")


# copy obj
sim_obj_i <- sim_obj_1

sim_obj_1 <- sim_obj_i
sim_obj_2 <- sim_obj_i
sim_obj_3 <- sim_obj_i
sim_obj_4 <- sim_obj_i
sim_obj_5 <- sim_obj_i

reset <- function() {
  sim_obj_1 <- sim_obj_i
  sim_obj_2 <- sim_obj_i
  sim_obj_3 <- sim_obj_i
  sim_obj_4 <- sim_obj_i
  sim_obj_5 <- sim_obj_i
  return('Sim_obj_1-5 have been reset')
}
reset()


#########################################
#########################################
#### finals without generating holes ####
#########################################
#########################################

#################
# max prob of 1 #
#################
sim_obj_1 <- GenerateCellPositivity(sim_obj_1, k = 10,
                                    sdmin = 0.5, sdmax = 2,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_2 <- GenerateCellPositivity(sim_obj_2, k = 10,
                                    sdmin = 0.5, sdmax = 1.5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_3 <- GenerateCellPositivity(sim_obj_3, k = 10,
                                    sdmin = 0.5, sdmax = 1,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_4 <- GenerateCellPositivity(sim_obj_4, k = 10,
                                    sdmin = 0.5, sdmax = 0.75,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_5 <- GenerateCellPositivity(sim_obj_5, k = 10,
                                    sdmin = 0.5, sdmax = 0.625,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)

PlotSimulation(sim_obj_1, which = 1, what = "whole core")
PlotSimulation(sim_obj_2, which = 1, what = "whole core")
PlotSimulation(sim_obj_3, which = 1, what = "whole core")
PlotSimulation(sim_obj_4, which = 1, what = "whole core")
PlotSimulation(sim_obj_5, which = 1, what = "whole core")


#############################
# recreating with 1 cluster #
#############################
sim_obj_1 <- GenerateCellPositivity(sim_obj_1, k = 1,
                                    sdmin = 0.5, sdmax = 2,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_2 <- GenerateCellPositivity(sim_obj_2, k = 1,
                                    sdmin = 0.5, sdmax = 1.5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_3 <- GenerateCellPositivity(sim_obj_3, k = 1,
                                    sdmin = 0.5, sdmax = 1,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_4 <- GenerateCellPositivity(sim_obj_4, k = 1,
                                    sdmin = 0.5, sdmax = 0.75,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_5 <- GenerateCellPositivity(sim_obj_5, k = 1,
                                    sdmin = 0.5, sdmax = 0.625,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)


PlotSimulation(sim_obj_1, which = 1, what = "whole core")
PlotSimulation(sim_obj_2, which = 1, what = "whole core")
PlotSimulation(sim_obj_3, which = 1, what = "whole core")
PlotSimulation(sim_obj_4, which = 1, what = "whole core")
PlotSimulation(sim_obj_5, which = 1, what = "whole core")


##################################
##################################
#### changing # of cell types ####
##################################
##################################
# create obj
sim_obj_1 <- CreateSimulationObject(sims = 1, cell_types = 1, window = custom_window)
sim_obj_2 <- CreateSimulationObject(sims = 1, cell_types = 2, window = custom_window)
sim_obj_3 <- CreateSimulationObject(sims = 1, cell_types = 3, window = custom_window)
sim_obj_4 <- CreateSimulationObject(sims = 1, cell_types = 4, window = custom_window)
sim_obj_5 <- CreateSimulationObject(sims = 1, cell_types = 5, window = custom_window)

# generate spatial pattern
sim_obj_1 <- GenerateSpatialPattern(sim_obj_1, 25)
sim_obj_2 <- GenerateSpatialPattern(sim_obj_2, 25)
sim_obj_3 <- GenerateSpatialPattern(sim_obj_3, 25)
sim_obj_4 <- GenerateSpatialPattern(sim_obj_4, 25)
sim_obj_5 <- GenerateSpatialPattern(sim_obj_5, 25)

plot(sim_obj_1, what = "Patterns", ncol = 1, nrow = 1, which = 1)
plot(sim_obj_2, what = "Patterns", ncol = 1, nrow = 1, which = 1)
plot(sim_obj_3, what = "Patterns", ncol = 1, nrow = 1, which = 1)
plot(sim_obj_4, what = "Patterns", ncol = 1, nrow = 1, which = 1)
plot(sim_obj_5, what = "Patterns", ncol = 1, nrow = 1, which = 1)

# generate tissue
sim_obj_1 <- GenerateTissue(sim_obj_1, density_heatmap = T, step_size = 0.1, cores = 1)
sim_obj_2 <- GenerateTissue(sim_obj_2, density_heatmap = T, step_size = 0.1, cores = 1)
sim_obj_3 <- GenerateTissue(sim_obj_3, density_heatmap = T, step_size = 0.1, cores = 1)
sim_obj_4 <- GenerateTissue(sim_obj_4, density_heatmap = T, step_size = 0.1, cores = 1)
sim_obj_5 <- GenerateTissue(sim_obj_5, density_heatmap = T, step_size = 0.1, cores = 1)


PlotSimulation(sim_obj_1, which = 1, ncol = 2, nrow = 2, what = "tissue heatmap")
PlotSimulation(sim_obj_2, which = 1, ncol = 2, nrow = 2, what = "tissue heatmap")
PlotSimulation(sim_obj_3, which = 1, ncol = 2, nrow = 2, what = "tissue heatmap")
PlotSimulation(sim_obj_4, which = 1, ncol = 2, nrow = 2, what = "tissue heatmap")
PlotSimulation(sim_obj_5, which = 1, ncol = 2, nrow = 2, what = "tissue heatmap")


# finals
sim_obj_1 <- GenerateCellPositivity(sim_obj_1, k = 1,
                                    sdmin = 0.5, sdmax = 2,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_2 <- GenerateCellPositivity(sim_obj_2, k = 1,
                                    sdmin = 0.5, sdmax = 2,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_3 <- GenerateCellPositivity(sim_obj_3, k = 1,
                                    sdmin = 0.5, sdmax = 2,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_4 <- GenerateCellPositivity(sim_obj_4, k = 1,
                                    sdmin = 0.5, sdmax = 2,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)
sim_obj_5 <- GenerateCellPositivity(sim_obj_5, k = 1,
                                    sdmin = 0.5, sdmax = 2,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                    shift = 1)


PlotSimulation(sim_obj_1, which = 1, what = "whole core")
PlotSimulation(sim_obj_2, which = 1, what = "whole core")
PlotSimulation(sim_obj_3, which = 1, what = "whole core")
PlotSimulation(sim_obj_4, which = 1, what = "whole core")
PlotSimulation(sim_obj_5, which = 1, what = "whole core")