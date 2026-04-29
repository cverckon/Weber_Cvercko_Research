library(scSpatialSIM)
library(ggplot2)
library(ggplotify)
library(patchwork)
library(tidyr)
set.seed(100)


maxs <- seq(9, 6.5, -0.5)
mins <- seq(1, 3.5, 0.5)
kdis <- maxs - mins

objs <- vector('list', 6)
pl <- vector('list', 6)


for (i in 1:length(objs)) {
  objs[[i]] <- vector('list', 20)
}

for (i in 1:length(maxs)) {
  for(j in 1:length(objs[[i]])) {
    custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
    # create obj
    objs[[i]][[j]] <- CreateSimulationObject(sims = 1, cell_types = 1, window = custom_window)
    
    # generate spatial pattern
    objs[[i]][[j]] <- GenerateSpatialPattern(objs[[i]][[j]], 25)
    
    # generate tissue
    objs[[i]][[j]] <- GenerateTissue(objs[[i]][[j]], xmin= mins[i], ymin= mins[i], xmax= maxs[i], ymax= maxs[i], step_size = 0.1, cores = 1)
    
    objs[[i]][[j]] <- GenerateCellPositivity(objs[[i]][[j]], k = 1, xmin= mins[i], ymin= mins[i], xmax= maxs[i], ymax= maxs[i],
                                             sdmin = 0.5, sdmax = 2,
                                             step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                             shift = 1)  
    pl[[i]][[j]] <- as.ggplot(PlotSimulation(objs[[i]][[j]], which = 1, what = "whole core"))
  }
}

##################
## OBJS CREATED ##
##################

saveRDS(objs, 'GenTissXY_OBJS.rds')
saveRDS(pl, 'GenTissXY_pl.rds')
objs <- readRDS('GenTissXY_OBJS.rds')
pl <- readRDS('GenTissXY_pl.rds')


maxs <- seq(9, 6.5, -0.5)
mins <- seq(1, 3.5, 0.5)
kdis <- maxs - mins


pf <- (pl[[1]][[1]] + pl[[2]][[1]]) / (pl[[3]][[1]] + pl[[4]][[1]]) / (pl[[5]][[1]] + pl[[6]][[1]])
ggsave("tissue_XY.png", plot = pf, dpi = 300, scale = 2)

################################
################################
####  CALCULATE FREQUENCIES ####
################################
################################

ca <- vector('list', length(objs))
cpos <- vector('list', length(objs))
cfreq <- vector('list', length(objs))

for (i in 1:length(objs)) {
  ca[[i]] <- vector('list', length(objs[[i]]))
  cfreq[[i]] <- vector('list', length(objs[[i]]))
  cpos[[i]] <- vector('list', length(objs[[i]]))
}


for (i in 1:length(objs)) { 
  for(j in 1:length(objs[[i]])) {
    ca[[i]][[j]] <- CreateSpatialList(objs[[i]][[j]])[[1]][3:4]
    cfreq[[i]][[j]] <- sum(ca[[i]][[j]][2] == 1 & ca[[i]][[j]][1] == 'Tissue 1' )/sum(ca[[i]][[j]][1] == 'Tissue 1')
    cpos[[i]][[j]] <- sum(ca[[i]][[j]][2] == 1 & ca[[i]][[j]][1] == 'Tissue 1' )
  }
}

f1 <- vector('list', length(objs))
p1 <- vector('list', length(objs))
v.f1 <- vector('list', length(objs))
v.p1 <- vector('list', length(objs))
for (i in 1:length(f1)) {
  f1[[i]] <- unlist(cfreq[[i]])
  p1[[i]] <- unlist(cpos[[i]])
  v.f1[[i]] <- var(f1[[i]])
  v.p1[[i]] <- var(p1[[i]])
}


# create dataframes for pos and freq
d.f1 <- as.data.frame(f1)
d.p1 <- as.data.frame(p1)

cn <- list(paste0('KernelDistances=', kdis))[[1]]
colnames(d.f1) <- cn
colnames(d.p1) <- cn

# convert to long frames

d.f1.long <- pivot_longer(d.f1,
                          cols = all_of(cn),
                          names_to = "kdis",
                          values_to = "Frequencies")

d.p1.long <- pivot_longer(d.p1,
                          cols = all_of(cn),
                          names_to = "kdis",
                          values_to = "Positives")

# create boxplots

freq.box <- ggplot(d.f1.long, aes(y= 0, x= Frequencies)) +
  geom_boxplot() + geom_jitter(height= 0.000001, size= 1.5, color= 'black', alpha= 0.4) +
  facet_wrap(~ kdis, ncol = 1) +
  theme_bw();freq.box
ggsave("GenTiss_xy_box_freq.png", plot = freq.box, dpi = 300, scale = 2)

pos.box <- ggplot(d.p1.long, aes(y= 0, x= Positives)) +
  geom_boxplot() + geom_jitter(height= 0.000001, size= 1.5, color= 'black', alpha= 0.4) +
  facet_wrap(~ kdis, ncol = 1) +
  theme_bw();pos.box
ggsave("GenTiss_xy_box_pos.png", plot = pos.box, dpi = 300, scale = 2)


##################################
##################################
####  MAXIMUM SD WERE VARIED  ####
##################################
##################################

param <- ExtractParameters(objs[[1]][[1]], 'All')
pv <- vector('list', length(param))
for (i in 1:length(param)){ 
  pv[[i]] <- unlist(param[[i]])
}
names(pv) <- names(param)
pv;pd

maxs
mins
dis



