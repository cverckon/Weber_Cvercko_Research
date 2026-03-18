#########################################
### varying sdmax in GenerateTissue() ###
#########################################
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


# copy obj
sim_obj_i <- sim_obj_1

sim_obj_1 <- sim_obj_i
sim_obj_2 <- sim_obj_i
sim_obj_3 <- sim_obj_i
sim_obj_4 <- sim_obj_i
sim_obj_5 <- sim_obj_i

objs <- c(sim_obj_1, sim_obj_2, sim_obj_3, sim_obj_4, sim_obj_5)

# generate tissue
maxs <- c(2, 1.5, 1, 0.75, 0.51)
for (i in 1:length(objs)) {
  objs[[i]] <- GenerateTissue(objs[[i]], density_heatmap = T, step_size = 0.1, cores = 1,
                                      sdmin = 0.5, sdmax = maxs[i])
}

sim_obj_1 <- GenerateTissue(sim_obj_1, density_heatmap = T, step_size = 0.1, cores = 1)
summary(sim_obj_1)

pl <- vector('list', 5)
for (i in 1:length(objs)) {
 pl[[i]] <- as.ggplot(PlotSimulation(objs[[i]], which = 1, ncol = 2, nrow = 2, what = "tissue heatmap")) 
}

#tmaps <- (pl[[1]] + pl[[2]] + plot_layout(widths = c(1, 1))) / 
  #(pl[[3]] + plot_layout(widths = c(1, 1))) / 
     #(pl[[4]] + pl[[5]] + plot_layout(widths = c(1, 1)))

tmaps <- pl[[1]] + pl[[2]] + pl[[3]] + pl[[4]] + pl[[5]] + plot_layout(ncol = 2)
ggsave("tissu-heatmaps.png", plot = tmaps, dpi = 300, scale = 2)

#############################
# recreating with 1 cluster #
#############################
for (i in 1:length(objs)) {
  objs[[i]] <- GenerateCellPositivity(objs[[i]], k = 1,
                                      sdmin = 0.5, sdmax = 2,
                                      density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                      shift = 1)
}

pl.2 <- vector('list', 5)
for (i in 1:length(pl.2)) {
  pl.2[[i]] <- as.ggplot(PlotSimulation(objs[[i]], which = 1, what = "whole core")) 
}
############################
##  calculate frequencies ##
############################

ca <- vector('list', length(objs))
cpos <- vector('list', length(objs))
cfreq <- vector('list', length(objs))
for (i in 1:length(objs)) {
  obj <- CreateSpatialList(objs[[i]])$`Spatial Data 1`
  index <- which(obj$`Tissue Assignment` == 'Tissue 1')
  ca[[i]] <- CreateSpatialList(objs[[i]])$`Spatial Data 1`$`Cell 1 Assignment`[index]
  cfreq[i] <- sum(ca[[i]] == 1)/length(ca[[i]])
  cpos[i] <- sum(ca[[i]] == 1)
}

ca2 <- vector('list', length(objs))
cpos2 <- vector('list', length(objs))
cfreq2 <- vector('list', length(objs))
for (i in 1:length(objs)) { 
  ca2[[i]] <- CreateSpatialList(objs[[i]])$`Spatial Data 1`$`Cell 1 Assignment`
  cfreq2[i] <- sum(ca2[[i]] == 1)/length(ca2[[i]])
  cpos2[i] <- sum(ca2[[i]] == 1)
}

f1 <- unlist(cfreq)
p1 <- unlist(cpos)
v.f1 <- var(f1)
v.p1 <- var(p1)
f2 <- unlist(cfreq2)
p2 <- unlist(cpos2)
v.f2 <- var(f2)
v.p2 <- var(p2)

v.f1; v.p1
v.f2; v.p2

pf <- (pl.2[[1]] + pl.2[[2]]) / pl.2[[3]] / (pl.2[[4]] + pl.2[[5]]); pf
ggsave("gen-tiss.png", plot = pf, dpi = 300, scale = 2)
