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


PlotSimulation(sim_obj_1, which = 1, ncol = 2, nrow = 2, what = "tissue heatmap")


# copy obj

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

objs <- c(sim_obj_1, sim_obj_2, sim_obj_3, 
           sim_obj_4, sim_obj_5)

ca <- vector('list', length(objs))
cpos <- vector('list', length(objs))
cfreq <- vector('list', length(objs))


for (i in 1:length(objs)) { 
  ca[[i]] <- CreateSpatialList(objs[[i]])$`Spatial Data 1`$`Cell 1 Assignment`
  cfreq[i] <- sum(ca[[i]] == 1)/length(ca[[i]])
  cpos[i] <- sum(ca[[i]] == 1)
}





for (i in 1:length(objs)){
  if(
    length(unique((ca[[i]]))) != 2
    ){print('error') 
  }
  }


counts <- data.frame(Positive = unlist(cpos), Frequencies = unlist(cfreq))
rownames(counts) <- list(paste0('Obj_', 1:length(objs)))[[1]]
counts











######### old

c1a <- slot(sim_obj_1, 'Spatial Files')[[1]]$`Cell 1 Assignment`
c1pos <- length(which(c1a == 'Positive'))
c1neg <- length(which(c1a != 'Positive'))
c1freq <- c1pos/c1neg

c2a <- slot(sim_obj_2, 'Spatial Files')[[1]]$`Cell 1 Assignment`
c2pos <- length(which(c2a == 'Positive'))
c2neg <- length(which(c2a != 'Positive'))
c2freq <- c2pos/c2neg

c3a <- slot(sim_obj_3, 'Spatial Files')[[1]]$`Cell 1 Assignment`
c3pos <- length(which(c3a == 'Positive'))
c3neg <- length(which(c3a != 'Positive'))
c3freq <- c3pos/c3neg

c4a <- slot(sim_obj_4, 'Spatial Files')[[1]]$`Cell 1 Assignment`
c4pos <- length(which(c4a == 'Positive'))
c4neg <- length(which(c4a != 'Positive'))
c4freq <- c4pos/c4neg

c5a <- slot(sim_obj_5, 'Spatial Files')[[1]]$`Cell 1 Assignment`
c5pos <- length(which(c5a == 'Positive'))
c5neg <- length(which(c5a != 'Positive'))
c5freq <- c5pos/c5neg

counts <- as.table(cbind(c(c1pos, c2pos, c3pos, c4pos, c5pos), 
                         c(c1neg, c2neg, c3neg, c4neg, c5neg), 
                         c(c1freq, c2freq, c3freq, c4freq, c5freq)
                         ))
colnames(counts) <-  c('positive', 'negative', 'freq')
rownames(counts) <- c('obj1', 'obj2', 'obj3', 'obj4', 'obj5')
counts

counts <- data.frame(cbind(cbind(cpos), cbind(cfreq)))
rownames(counts) <- list(paste0('Obj_', 1:length(objs)))[[1]]
colnames(counts) <- c('Positive', 'Frequencies')
counts

################################

CreateSpatialList(single_df = FALSE)
spat_data_distribution = GenerateDistributions(sim_obj_1)






############
##  Final ##
############
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


PlotSimulation(sim_obj_1, which = 1, ncol = 2, nrow = 2, what = "tissue heatmap")

sim_obj_i <- sim_obj_1

sim_obj_1 <- sim_obj_i
sim_obj_2 <- sim_obj_i
sim_obj_3 <- sim_obj_i
sim_obj_4 <- sim_obj_i
sim_obj_5 <- sim_obj_i

objs <- c(sim_obj_1, sim_obj_2, sim_obj_3, sim_obj_4, sim_obj_5)

maxs <- c(2, 1.5, 1, 0.76, 0.51)
for (i in 1:length(objs)) {
  objs[[i]] <- GenerateCellPositivity(objs[[i]], k = 1,
                                      sdmin = 0.5, sdmax = maxs[i],
                                      density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                      shift = 1)
}





pl <- vector('list', length(objs))
for (i in 1: length(objs)) {
  pl[[i]] <- PlotSimulation(objs[[i]], which = 1, what = "whole core")
}

pl[[5]]; pl[[4]]; pl[[3]]; pl[[2]]; pl[[1]]



ca <- vector('list', length(objs))
cpos <- vector('list', length(objs))
cfreq <- vector('list', length(objs))


for (i in 1:length(objs)) { 
  ca[[i]] <- CreateSpatialList(objs[[i]])$`Spatial Data 1`$`Cell 1 Assignment`
  cfreq[i] <- sum(ca[[i]] == 1)/length(ca[[i]])
  cpos[i] <- sum(ca[[i]] == 1)
}

counts <- data.frame(Positive = unlist(cpos), Frequencies = unlist(cfreq))
rownames(counts) <- list(paste0('Obj_', 1:length(objs)))[[1]]
counts
