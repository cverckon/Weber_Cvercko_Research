library(scSpatialSIM)
library(ggplot2)
library(ggplotify)
library(patchwork)
set.seed(800)

objs <- vector('list', 6)
pl <- vector('list', 6)

for (i in 1:length(objs)) {
  objs[[i]] <- vector('list', 20)
}

for (i in 1:length(s)) {
  for(j in 1:length(objs[[i]])) {
    custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
    # create obj
    objs[[i]][[j]] <- CreateSimulationObject(sims = 1, cell_types = 1, window = custom_window)
    
    # generate spatial pattern
    objs[[i]][[j]] <- GenerateSpatialPattern(objs[[i]][[j]], 25)
    
    # generate tissue
    objs[[i]][[j]] <- GenerateTissue(objs[[i]][[j]], step_size = 0.1, cores = 1)
    
    objs[[i]][[j]] <- GenerateCellPositivity(objs[[i]][[j]], k = 1, xmin= 2.5, ymin= 2.5, xmax= 7.5, ymax= 7.5,
                                             sdmin = 0.5, sdmax = 2,
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


d.f1 <- as.data.frame(f1)
d.p1 <- as.data.frame(p1)

c <- list(paste0('SimSet_', 1:length(objs)))[[1]]
colnames(d.f1) <- c
colnames(d.p1) <- c

fb1 <- ggplot(d.f1, aes(x= SimSet_1, y= 0)) + geom_boxplot() + xlim(0,0.1) + geom_jitter(height= 0.0001)
fb2 <- ggplot(d.f1, aes(x= SimSet_2, y= 0)) + geom_boxplot() + xlim(0,0.1) + geom_jitter(height= 0.0001)
fb3 <- ggplot(d.f1, aes(x= SimSet_3, y= 0)) + geom_boxplot() + xlim(0,0.1) + geom_jitter(height= 0.0001)
fb4 <- ggplot(d.f1, aes(x= SimSet_4, y= 0)) + geom_boxplot() + xlim(0,0.1) + geom_jitter(height= 0.0001)
fb5 <- ggplot(d.f1, aes(x= SimSet_5, y= 0)) + geom_boxplot() + xlim(0,0.1) + geom_jitter(height= 0.0001)
fb6 <- ggplot(d.f1, aes(x= SimSet_6, y= 0)) + geom_boxplot() + xlim(0,0.1) + geom_jitter(height= 0.0001)

freq.box <- fb1 / fb2 / fb3 / fb4 / fb5 / fb6
ggsave("seeds_control_freq_box.png", plot = freq.box, dpi = 300, scale = 2)


pb1 <- ggplot(d.p1, aes(x= SimSet_1, y= 0)) + geom_boxplot() + geom_jitter(height= 0.0001) + xlim(0, 200)
pb2 <- ggplot(d.p1, aes(x= SimSet_2, y= 0)) + geom_boxplot()  + geom_jitter(height= 0.0001) + xlim(0, 200)
pb3 <- ggplot(d.p1, aes(x= SimSet_3, y= 0)) + geom_boxplot() + geom_jitter(height= 0.0001) + xlim(0, 200)
pb4 <- ggplot(d.p1, aes(x= SimSet_4, y= 0)) + geom_boxplot() + geom_jitter(height= 0.0001) + xlim(0, 200)
pb5 <- ggplot(d.p1, aes(x= SimSet_5, y= 0)) + geom_boxplot() + geom_jitter(height= 0.0001) + xlim(0, 200)
pb6 <- ggplot(d.p1, aes(x= SimSet_6, y= 0 )) + geom_boxplot() + geom_jitter(height= 0.0001) + xlim(0, 200)

pos.box <- pb1 / pb2 / pb3 / pb4 / pb5 / pb6
ggsave("seeds_control_pos_box.png", plot = pos.box, dpi = 300, scale = 2)


###################################################################
###################################################################
####  The only Difference accross trials were the seed values  ####
###################################################################
###################################################################

param <- ExtractParameters(objs[[1]], 'All')
pv <- vector('list', length(param))
for (i in 1:length(param)){ 
  pv[[i]] <- unlist(param[[i]])
}
names(pv) <- names(param)
pv
s


