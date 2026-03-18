library(scSpatialSIM)
library(ggplot2)
library(ggplotify)
library(patchwork)
set.seed(100) #reproducibility

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


PlotSimulation(sim_obj_1, which = 1, ncol = 2, nrow = 2, what = "tissue heatmap")


# copy obj


sim_obj_i <- sim_obj_1

sim_obj_1 <- sim_obj_i
sim_obj_2 <- sim_obj_i
sim_obj_3 <- sim_obj_i
sim_obj_4 <- sim_obj_i
sim_obj_5 <- sim_obj_i

objs <- c(sim_obj_1, sim_obj_2, sim_obj_3, sim_obj_4, sim_obj_5)


#############################
# recreating with 1 cluster #
#############################


maxs <- c(2, 1.5, 1, 0.76, 0.51)
for (i in 1:length(objs)) {
  objs[[i]] <- GenerateCellPositivity(objs[[i]], k = 1,
                                      sdmin = 0.5, sdmax = maxs[i],
                                      density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                      shift = 1)
}

pl1 <- as.ggplot(PlotSimulation(objs[[1]], which = 1, what = "whole core"))
pl2 <- as.ggplot(PlotSimulation(objs[[2]], which = 1, what = "whole core"))
pl3 <- as.ggplot(PlotSimulation(objs[[3]], which = 1, what = "whole core"))
pl4 <- as.ggplot(PlotSimulation(objs[[4]], which = 1, what = "whole core"))
pl5 <- as.ggplot(PlotSimulation(objs[[5]], which = 1, what = "whole core"))


sim_obj_1@Cells

slotNames(sim_obj_1)
slot(sim_obj_1, 'Cells')

CreateSpatialList(objs[[1]])
############################
##  calculate frequencies ##
############################

ca <- vector('list', length(objs))
cpos <- vector('list', length(objs))
cfreq <- vector('list', length(objs))


for (i in 1:length(objs)) { 
  ca[[i]] <- CreateSpatialList(objs[[i]])$`Spatial Data 1`$`Cell 1 Assignment`
  cfreq[i] <- sum(ca[[i]] == 1)/length(ca[[i]])
  cpos[i] <- sum(ca[[i]] == 1)
}


f1 <- unlist(cfreq);f1
p1 <- unlist(cpos);p1

v.f1 <- var(f1); v.f1
v.p1 <- var(p1); v.p1



pf <- (pl1 + pl2) / pl3 / (pl4 + pl5); pf

ggsave("genpos-plots.png", plot = pf, dpi = 300, scale = 2)
