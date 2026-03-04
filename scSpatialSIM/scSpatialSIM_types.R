library(scSpatialSIM)
set.seed(333) #reproducibility

custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

# create obj

objs <- vector('list', 5)
for (i in 1:length(objs)) {
  objs[[i]] <- CreateSimulationObject(sims = 1, cell_types = i, window = custom_window)
}


# generate spatial pattern
for (i in 1:length(objs)) {
  objs[[i]] <- GenerateSpatialPattern(objs[[i]], 50)
}


summary(objs[[1]])
plot(objs[[1]], what = "Patterns", ncol = 1, nrow = 1, which = 1) 


# generate tissue

for (i in 1:length(objs)) {
  objs[[i]] <- GenerateTissue(objs[[i]], density_heatmap = T, step_size = 0.1, cores = 1)
}


summary(objs[[1]])
PlotSimulation(objs[[1]], which = 1, ncol = 2, nrow = 2, what = "tissue heatmap")


#################################
# Create Probability dataframes #
# as number of cell type        #
# increases the probability of  #
# assignment will decrease.     #
#################################

p <- vector('list', length(objs))
for (i in 1:length(p)) {
  tp.max <- vector('list', i)
  for (i in 1:length(tp.max)) {
    min <- rep(0, i)
    tp.max[[i]] <- 1/i 
  }
  
  p[[i]] <- data.frame(min, max = unlist(tp.max))
}




for (i in 1:length(objs)) {
  objs[[i]] <- GenerateCellPositivity(objs[[i]], k = 1, xmin = 2.5, xmax = 7.5, ymin = 2.5, ymax = 7.5,
                                      sdmin = 0.5, sdmax = 2,
                                      density_heatmap = T, step_size = 0.1, cores = 1, probs = p[[i]],
                                      shift = 1)
}

##########################
# no_kernel = TRUE trial #
##########################

for (i in 1:length(objs)) {
  objs[[i]] <- GenerateCellPositivity(objs[[i]], no_kernel = TRUE,
                                      density_heatmap = T, step_size = 0.1, cores = 1, probs = p[[i]],
                                      shift = 1)
}

#####################################################################################################################################
pl <- vector('list', length(objs))
for (i in 1: length(objs)) {
  pl[[i]] <- PlotSimulation(objs[[i]], which = 1, what = "whole core")
}

pl[[5]]; pl[[4]]; pl[[3]]; pl[[2]]; pl[[1]]


cfs <- vector('list', length(objs))
cps <- vector('list', length(objs))
for (i in 1:length(objs)) { 
  ct <- CreateSpatialList(objs[[i]])$`Spatial Data 1`
  
  ca <- vector('list', ncol(ct) - 3)
  cfreq <- vector('list', length(ca))
  cpos <- vector('list', length(ca))
  for (i in 1:length(ca)) {
    ca[[i]] <- ct[i+3]
    cfreq[[i]] <- sum(ca[[i]] == 1)/nrow(ca[[i]])
    cpos[[i]] <- sum(ca[[i]] == 1)
  }
  cfs[[i]] <- cfreq
  cps[[i]] <- cpos
  }

counts <- data.frame(Positive = unlist(cpos), Frequencies = unlist(cfreq))
rownames(counts) <- list(paste0('Obj_', 1:length(objs)))[[1]]
counts