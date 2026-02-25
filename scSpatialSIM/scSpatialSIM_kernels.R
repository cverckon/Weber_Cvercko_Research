library(scSpatialSIM)
set.seed(333) #reproducibility

custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

# create obj
sim_obj_1 <- CreateSimulationObject(sims = 1, cell_types = 1, window = custom_window) %>%
  GenerateSpatialPattern() %>%
  GenerateTissue(density_heatmap = TRUE, step_size = 0.1, cores = 1)
summary(sim_obj_1)


library(scSpatialSIM)
set.seed(333) #reproducibility

custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

# create obj
sim_obj_1 <- CreateSimulationObject(sims = 1, cell_types = 1, window = custom_window)
summary(sim_obj_1)


# generate spatial pattern
sim_obj_1 <- GenerateSpatialPattern(sim_obj_1, 25)
summary(sim_obj_1)

plot(sim_obj_1, what = "Patterns", ncol = 1, nrow = 1, which = 1) 


# generate tissue
sim_obj_1 <- GenerateTissue(sim_obj_1, density_heatmap = T, step_size = 0.1, cores = 1)
summary(sim_obj_1)



for (i in 1:length(objs)) {
  objs[[i]] <- sim_obj_i
}

sim_obj_i <- sim_obj_1

sim_obj_1 <- sim_obj_i
sim_obj_2 <- sim_obj_i
sim_obj_3 <- sim_obj_i
sim_obj_4 <- sim_obj_i
sim_obj_5 <- sim_obj_i

objs <- c(sim_obj_1, sim_obj_2, sim_obj_3, sim_obj_4, sim_obj_5)



##################
## change sdmin ##
##################

mins <- c(0.5, 0.75, 1, 1.5, 1.9)
for (i in 1:length(objs)) {
  objs[[i]] <- GenerateCellPositivity(objs[[i]], k = 5,
                         sdmin = mins[i], sdmax = 2,
                         density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                         shift = 1)
}

PlotSimulation(objs[[1]], which = 1, what = 'whole core')
PlotSimulation(objs[[2]], which = 1, what = 'whole core')
PlotSimulation(objs[[3]], which = 1, what = 'whole core')
PlotSimulation(objs[[4]], which = 1, what = 'whole core')
PlotSimulation(objs[[5]], which = 1, what = 'whole core')

for (i in 1:length(objs)) {
  PlotSimulation(objs[[i]], which = 1, what = 'whole core')
}


mins <- c(0.5, 0.75, 1, 1.5, 1.9)
for (i in 1:length(objs)) {
  objs[[i]] <- GenerateCellPositivity(objs[[i]], k = 1,
                                      sdmin = mins[i], sdmax = 2, ymin = 5, ymax = 5, xmin = 5, xmax = 5,
                                      density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                      shift = 1)
}

PlotSimulation(objs[[1]], which = 1, what = 'whole core')
PlotSimulation(objs[[2]], which = 1, what = 'whole core')
PlotSimulation(objs[[3]], which = 1, what = 'whole core')
PlotSimulation(objs[[4]], which = 1, what = 'whole core')
PlotSimulation(objs[[5]], which = 1, what = 'whole core')


