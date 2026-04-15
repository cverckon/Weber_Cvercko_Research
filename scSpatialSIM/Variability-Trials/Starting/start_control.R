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

for (i in 1:length(objs)) {
  for(j in 1:length(objs[[i]])) {
    custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
    # create obj
    objs[[i]][[j]] <- CreateSimulationObject(sims = 1, cell_types = 1, window = custom_window)
    
    # generate spatial pattern
    objs[[i]][[j]] <- GenerateSpatialPattern(objs[[i]][[j]], 25)
    
    # generate tissue
    objs[[i]][[j]] <- GenerateTissue(objs[[i]][[j]], step_size = 0.1, cores = 1)
    
    objs[[i]][[j]] <- GenerateCellPositivity(objs[[i]][[j]], k = 1,
                                             sdmin = 0.5, sdmax = 2,
                                             step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                             shift = 1)  
    pl[[i]][[j]] <- as.ggplot(PlotSimulation(objs[[i]][[j]], which = 1, what = "whole core"))
  }
}


##################
## OBJS CREATED ##
##################

saveRDS(objs, 'StartControlOBJS.rds')
objs <- readRDS('StartControlOBJS.rds')


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


# create dataframes for pos and freq
d.f1 <- as.data.frame(f1)
d.p1 <- as.data.frame(p1)



cn <- list(paste0('sim_set', 1:6))[[1]]
colnames(d.f1) <- cn
colnames(d.p1) <- cn

# convert to long frames

d.f1.long <- pivot_longer(d.f1,
                          cols = all_of(cn),
                          names_to = "sim_set",
                          values_to = "Frequencies")

d.p1.long <- pivot_longer(d.p1,
                          cols = all_of(cn),
                          names_to = "sim_set",
                          values_to = "Frequencies")

# create boxplots

freq.box <- ggplot(d.f1.long, aes(y= 0, x= Frequencies)) +
  geom_boxplot() + geom_jitter(height= 0.000001, size= 1.5, color= 'black', alpha= 0.4) +
  facet_wrap(~ sim_set, ncol = 1) +
  theme_bw();freq.box
ggsave("start_control_box_freq.png", plot = freq.box, dpi = 300, scale = 2)

pos.box <- ggplot(d.p1.long, aes(y= 0, x= Frequencies)) +
  geom_boxplot() + geom_jitter(height= 0.000001, size= 1.5, color= 'black', alpha= 0.4) +
  facet_wrap(~ sim_set, ncol = 1) +
  theme_bw();pos.box
ggsave("start_control_box_pos.png", plot = pos.box, dpi = 300, scale = 2)


###################################################
###################################################
####  There were no differences across trials  ####
###################################################
###################################################

param <- ExtractParameters(objs[[1]], 'All')
pv <- vector('list', length(param))
for (i in 1:length(param)){ 
  pv[[i]] <- unlist(param[[i]])
}
names(pv) <- names(param)
pv 



