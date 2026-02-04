#if (!require("devtools", quietly = TRUE))
  #install.packages("devtools")

#devtools::install_github("FridleyLab/mIFsim@v0.1.3.4")

library(scSpatialSIM)
set.seed(333) #reproducibility

custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

sim_obj_1 <- CreateSimulationObject(sims = 9, cell_types = 1, window = custom_window)
sim_obj_2 <- CreateSimulationObject(sims = 9, cell_types = 1, window = custom_window)


summary(sim_obj_1)
summary(sim_obj_2)

sim_obj_1 <- GenerateSpatialPattern(sim_obj_1, 100)
sim_obj_2 <- GenerateSpatialPattern(sim_obj_2, 1)
summary(sim_obj_1)
summary(sim_obj_2)



plot(sim_obj_1, what = "Patterns", ncol = 1, nrow = 1, which = 1)#print only first point pattern
plot(sim_obj_2, what = "Patterns", ncol = 1, nrow = 1, which = 1)


sim_obj_1 <- GenerateTissue(sim_obj_1, density_heatmap = T, step_size = 0.1, cores = 1)
sim_obj_2 <- GenerateTissue(sim_obj_2, density_heatmap = T, step_size = 0.1, cores = 1)
summary(sim_obj_1)
summary(sim_obj_2)


PlotSimulation(sim_obj_1, which = 1:4, ncol = 2, nrow = 2, what = "tissue heatmap")
PlotSimulation(sim_obj_2, which = 1:4, ncol = 2, nrow = 2, what = "tissue heatmap")



sim_obj_1 <- GenerateHoles(sim_obj_1, density_heatmap = T, step_size = 0.1, cores = 1)
sim_obj_2 <- GenerateHoles(sim_obj_2, density_heatmap = T, step_size = 0.1, cores = 1)
summary(sim_obj_1)
summary(sim_obj_2)

PlotSimulation(sim_obj_1, which = 1:8, ncol = 2, nrow = 2, what = "hole heatmap")
PlotSimulation(sim_obj_2, which = 1:8, ncol = 2, nrow = 2, what = "hole heatmap")


sim_obj_1 <- GenerateCellPositivity(sim_obj_1, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)
sim_obj_2 <- GenerateCellPositivity(sim_obj_2, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)
summary(sim_obj_1)
summary(sim_obj_2)

PlotSimulation(sim_obj_1, which = 1, what = "whole core")
PlotSimulation(sim_obj_2, which = 1, what = "whole core")


###################################
## Testing Combining of sim_objs ##
###################################


t1 <- sim_obj_1@Patterns[[1]][1:4, ]
t2 <- sim_obj_2@Patterns[[1]][1:4, ]
tbind <- rbind(sim_obj_1@Patterns[[1]][1:4, ], sim_obj_2@Patterns[[1]][1:4, ]);tbind
fbind <- rbind(sim_obj_1@Patterns[[1]][ , 1], sim_obj_2@Patterns[[1]][ , 1])


x <- rbind(sim_obj_1@Patterns[[1]], sim_obj_2@Patterns[[1]])





class(sim_obj_1@Patterns)
sim_obj_3@Patterns <- c(sim_obj_3@Patterns, sim_obj_2@Patterns)

length(sim_obj_1@Patterns[[1]][ ,1]) + length(sim_obj_2@Patterns[[1]][ ,1])
length(x[ , 1])



###############
## COMBINING ##
###############

sim_obj_3 <- sim_obj_1
for (i in 1:length(sim_obj_3@Patterns)) {
  sim_obj_3@Patterns[[i]] <- rbind(sim_obj_3@Patterns[[i]], sim_obj_2@Patterns[[i]])
}

length(sim_obj_3@Patterns[[1]][ ,1])
length(sim_obj_1@Patterns[[1]][ ,1]) + length(sim_obj_2@Patterns[[1]][ ,1])

PlotSimulation(sim_obj_3, which = 1, what = "whole core")

####################################################
## OR COMBINE EARLIER AND GENERATE THE REST AFTER ##
####################################################

library(scSpatialSIM)
set.seed(333) #reproducibility

custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

sim_obj_1 <- CreateSimulationObject(sims = 9, cell_types = 1, window = custom_window)
sim_obj_2 <- CreateSimulationObject(sims = 9, cell_types = 1, window = custom_window)
summary(sim_obj_1)
summary(sim_obj_2)


sim_obj_1 <- GenerateSpatialPattern(sim_obj_1, 100)
sim_obj_2 <- GenerateSpatialPattern(sim_obj_2, 1)
summary(sim_obj_1)
summary(sim_obj_2)
plot(sim_obj_1, what = "Patterns", ncol = 1, nrow = 1, which = 1)#print only first point pattern
plot(sim_obj_2, what = "Patterns", ncol = 1, nrow = 1, which = 1)

sim_obj_3 <- sim_obj_1
for (i in 1:length(sim_obj_3@Patterns)) {
  sim_obj_3@Patterns[[i]] <- rbind(sim_obj_3@Patterns[[i]], sim_obj_2@Patterns[[i]])
}

length(sim_obj_3@Patterns[[1]][ ,1])
length(sim_obj_1@Patterns[[1]][ ,1]) + length(sim_obj_2@Patterns[[1]][ ,1])

sim_obj_3 <- GenerateTissue(sim_obj_3, density_heatmap = T, step_size = 0.1, cores = 1)
summary(sim_obj_3)
PlotSimulation(sim_obj_3, which = 1:4, ncol = 2, nrow = 2, what = "tissue heatmap")

sim_obj_3 <- GenerateHoles(sim_obj_3, density_heatmap = T, step_size = 0.1, cores = 1)
summary(sim_obj_3)
PlotSimulation(sim_obj_3, which = 1:8, ncol = 2, nrow = 2, what = "hole heatmap")


sim_obj_3 <- GenerateCellPositivity(sim_obj_3, k = 4,
                                    sdmin = 3, sdmax = 5,
                                    density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                    shift = 1)
summary(sim_obj_3)
PlotSimulation(sim_obj_3, which = 1, what = "whole core")

##########################
## Bivariate Simulation ##
##########################

#set seed
set.seed(333)
#create the new object
bivariate_sim <- CreateSimulationObject(sims = 5, cell_types = 2) %>%
  #produce the point pattern
  GenerateSpatialPattern() %>%
  #make tissues
  GenerateTissue(density_heatmap = T, step_size = 0.1, cores = 1)

bivariate_sim_tmp  <- GenerateCellPositivity(bivariate_sim, k = 4,
                                            sdmin = 1, sdmax = 2,
                                            density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                            shift = 0)

PlotSimulation(bivariate_sim_tmp, which = 1, what = "whole core")

bivariate_sim_tmp  <- GenerateCellPositivity(bivariate_sim, k = 4,
                                            sdmin = 3, sdmax = 5,
                                            density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 0.1),
                                            shift = 1)

PlotSimulation(bivariate_sim_tmp, which = 1, what = "whole core")

bivariate_sim_tmp@Cells
class(bivariate_sim_tmp@Cells)
mode(bivariate_sim_tmp@Cells)

bivariate_sim_tmp@Patterns

bivariate_sim_tmp@`Spatial Files`[[1]][1:20, ]
