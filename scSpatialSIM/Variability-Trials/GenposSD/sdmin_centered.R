library(scSpatialSIM)
library(ggplot2)
library(ggplotify)
library(patchwork)
set.seed(100)

maxs <- c(1.5, 1.25, 1, 0.75, 62.5 , 0.5)
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
    objs[[i]][[j]] <- GenerateTissue(objs[[i]][[j]], step_size = 0.1, cores = 1)
    
    objs[[i]][[j]] <- GenerateCellPositivity(objs[[i]][[j]], k = 1, xmin= 2.5, ymin= 2.5, xmax= 7.5, ymax= 7.5,
                                             sdmin = maxs[i], sdmax = 2,
                                             step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                             shift = 1)  
    pl[[i]][[j]] <- as.ggplot(PlotSimulation(objs[[i]][[j]], which = 1, what = "whole core"))
  }
}




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
    cfreq[[i]][[j]] <- sum(ca[[i]][[j]][2] == 1 & ca[[i]][[j]][1] == 'Tissue 2' )/sum(ca[[i]][[j]][1] == 'Tissue 2')
    cpos[[i]][[j]] <- sum(ca[[i]][[j]][2] == 1 & ca[[i]][[j]][1] == 'Tissue 2' )
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

####### box plots
d.f1 <- as.data.frame(f1)
d.p1 <- as.data.frame(p1)

c <- list(paste0('S', s))[[1]]
colnames(d.f1) <- c
colnames(d.p1) <- c

fb1 <- ggplot(d.f1, aes(x= S100, y= 0)) + geom_boxplot() + xlim(0,0.16) + geom_jitter(height= 0.0001)
fb2 <- ggplot(d.f1, aes(x= S200, y= 0)) + geom_boxplot() + xlim(0,0.16) + geom_jitter(height= 0.0001)
fb3 <- ggplot(d.f1, aes(x= S300, y= 0)) + geom_boxplot() + xlim(0,0.16) + geom_jitter(height= 0.0001)
fb4 <- ggplot(d.f1, aes(x= S400, y= 0)) + geom_boxplot() + xlim(0,0.16) + geom_jitter(height= 0.0001)
fb5 <- ggplot(d.f1, aes(x= S500, y= 0)) + geom_boxplot() + xlim(0,0.16) + geom_jitter(height= 0.0001)
fb6 <- ggplot(d.f1, aes(x= S600, y= 0 )) + geom_boxplot() + xlim(0,0.16) + geom_jitter(height= 0.0001)

freq.box <- fb1 / fb2 / fb3 / fb4 / fb5 / fb6
ggsave("sdmin_box_freq.png", plot = freq.box, dpi = 300, scale = 2)

pb1 <- ggplot(d.p1, aes(x= SD1.5, y= 0)) + geom_boxplot() + geom_jitter(height= 0.0001) + xlim(0, 300)
pb2 <- ggplot(d.p1, aes(x= SD1.25, y= 0)) + geom_boxplot()  + geom_jitter(height= 0.0001) + xlim(0, 300)
pb3 <- ggplot(d.p1, aes(x= SD1, y= 0)) + geom_boxplot() + geom_jitter(height= 0.0001) + xlim(0, 300)
pb4 <- ggplot(d.p1, aes(x= SD0.75, y= 0)) + geom_boxplot() + geom_jitter(height= 0.0001) + xlim(0, 300)
pb5 <- ggplot(d.p1, aes(x= SD0.625, y= 0)) + geom_boxplot() + geom_jitter(height= 0.0001) + xlim(0, 300)
pb6 <- ggplot(d.p1, aes(x= SD0.5, y= 0 )) + geom_boxplot() + geom_jitter(height= 0.0001) + xlim(0, 300)

pos.box <- pb1 / pb2 / pb3 / pb4 / pb5 / pb6
ggsave("sdmin_box_pos.png", plot = pos.box, dpi = 300, scale = 2)

##################################
##################################
####  MINIMUM SD WERE VARIED  ####
##################################
##################################

param <- ExtractParameters(objs[[1]], 'All')
pv <- vector('list', length(param))
for (i in 1:length(param)){ 
  pv[[i]] <- unlist(param[[i]])
}
names(pv) <- names(param)
pv;pd

maxs

