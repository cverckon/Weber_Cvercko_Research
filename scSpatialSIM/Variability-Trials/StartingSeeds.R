############################
############################
############################
######                ######
######   Seed = 100   ######
######                ######
############################
############################
############################

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

PlotSimulation(objs[[1]], which = 1, what = "whole core")
PlotSimulation(objs[[2]], which = 1, what = "whole core")
PlotSimulation(objs[[3]], which = 1, what = "whole core")
PlotSimulation(objs[[4]], which = 1, what = "whole core")
PlotSimulation(objs[[5]], which = 1, what = "whole core")


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


f1 <- unlist(cfreq)
p1 <- unlist(cpos)

############################
############################
############################
######                ######
######   Seed = 200   ######
######                ######
############################
############################
############################

set.seed(200) #reproducibility

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

PlotSimulation(objs[[1]], which = 1, what = "whole core")
PlotSimulation(objs[[2]], which = 1, what = "whole core")
PlotSimulation(objs[[3]], which = 1, what = "whole core")
PlotSimulation(objs[[4]], which = 1, what = "whole core")
PlotSimulation(objs[[5]], which = 1, what = "whole core")


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

f2 <- unlist(cfreq)
p2 <- unlist(cpos)

############################
############################
############################
######                ######
######   Seed = 300   ######
######                ######
############################
############################
############################

set.seed(300) #reproducibility

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

PlotSimulation(objs[[1]], which = 1, what = "whole core")
PlotSimulation(objs[[2]], which = 1, what = "whole core")
PlotSimulation(objs[[3]], which = 1, what = "whole core")
PlotSimulation(objs[[4]], which = 1, what = "whole core")
PlotSimulation(objs[[5]], which = 1, what = "whole core")


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

f3 <- unlist(cfreq)
p3 <- unlist(cpos)



############################
############################
############################
######                ######
######   Seed = 400   ######
######                ######
############################
############################
############################

set.seed(400) #reproducibility

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

PlotSimulation(objs[[1]], which = 1, what = "whole core")
PlotSimulation(objs[[2]], which = 1, what = "whole core")
PlotSimulation(objs[[3]], which = 1, what = "whole core")
PlotSimulation(objs[[4]], which = 1, what = "whole core")
PlotSimulation(objs[[5]], which = 1, what = "whole core")


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

f4 <- unlist(cfreq)
p4 <- unlist(cpos)

############################
############################
############################
######                ######
######   Seed = 500   ######
######                ######
############################
############################
############################

set.seed(500) #reproducibility

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

PlotSimulation(objs[[1]], which = 1, what = "whole core")
PlotSimulation(objs[[2]], which = 1, what = "whole core")
PlotSimulation(objs[[3]], which = 1, what = "whole core")
PlotSimulation(objs[[4]], which = 1, what = "whole core")
PlotSimulation(objs[[5]], which = 1, what = "whole core")


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

f5 <- unlist(cfreq)
p5 <- unlist(cpos)

############################################

fs <- vector('list', 5)
fs[[1]] <- f1; fs[[2]] <- f2; fs[[3]] <- f3; fs[[4]] <- f4; fs[[5]] <- f5
ps <- vector('list', 5)
ps[[1]] <- p1; ps[[2]] <- p2; ps[[3]] <- p3; ps[[4]] <- p4; ps[[5]] <- p5

fs.avg <- vector('list', 5)
ps.avg <- vector('list', 5)
for (i in 1:length(fs)) {
  fs.avg[[i]] <- mean(fs[[i]])
  ps.avg[[i]] <- mean(ps[[i]])
  }

var.fs.avg <- var(unlist(fs.avg));var.fs.avg
var.ps.avg <- var(unlist(ps.avg));var.ps.avg


fs.obj <- vector('list', 5)
ps.obj <- vector('list', 5)
var.fs.obj <- vector('list', 5)
var.ps.obj <- vector('list', 5)
for (i in 1:length(fs.obj)) {
  fs.obj[[i]] <- vector('list', 5)
  ps.obj[[i]] <- vector('list', 5)
  for (j in 1:length(fs.obj[[i]])) {
    fs.obj[[i]][j] <- fs[[j]][i]
    ps.obj[[i]][j] <- ps[[j]][i]  
    }
  fs.obj[[i]] <- unlist(fs.obj[[i]])
  ps.obj[[i]] <- unlist(ps.obj[[i]])
  var.fs.obj[[i]] <- var(fs.obj[[i]])
  var.ps.obj[[i]] <- var(ps.obj[[i]])
}

v <- data.frame(sdMax = maxs, var_freq = unlist(var.fs.obj), var_pos = unlist(var.ps.obj))

counts <- data.frame(Positive = unlist(cpos), Frequencies = unlist(cfreq))
rownames(counts) <- list(paste0('Obj_', 1:length(objs)))[[1]]
counts



