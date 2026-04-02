library(scSpatialSIM)
library(ggplot2)
library(ggplotify)
library(patchwork)
set.seed(100)
ns <- 6

custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
# create obj
obj1 <- CreateSimulationObject(sims= ns, cell_types = 1, window = custom_window)
  
# generate spatial pattern
obj1 <- GenerateSpatialPattern(obj1, 25)
  
# generate tissue
obj1 <- GenerateTissue(obj1, density_heatmap = T, step_size = 0.1, cores = 1)
  
obj1 <- GenerateCellPositivity(obj1, k = 1, xmin= 2.5, ymin= 2.5, xmax= 7.5, ymax= 7.5,
                                      sdmin = 0.5, sdmax = 2,
                                      density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                      shift = 1)  

pl <- vector('list', ns)
for (i in 1:ns){
  pl[[i]] <- as.ggplot(PlotSimulation(obj1, which = i, what = "whole core"))
}

################################
################################
####  CALCULATE FREQUENCIES ####
################################
################################

  


ca <- vector('list', ns)
cpos <- vector('list', ns)
cfreq <- vector('list', ns)
for (i in 1:ns) { 
  ca[[i]] <- CreateSpatialList(obj1)[[i]]$`Cell 1 Assignment`
  cfreq[i] <- sum(ca[[i]] == 1)/length(ca[[i]])
  cpos[i] <- sum(ca[[i]] == 1)
}

f1 <- unlist(cfreq)
p1 <- unlist(cpos)
v.f1 <- var(f1)
v.p1 <- var(p1)

f1 <- c(f1, v.f1)
p1 <- c(p1, v.p1)
rs <- list(paste0('Simulation ', 1:ns))[[1]]
rs <- c(rs, 'VARIANCE')

counts <- data.frame(Positive = p1, Simulations= rs, Frequencies = f1)
counts

counts.2 <- counts[-(nrow(counts)), ]
freq.bars <- ggplot(counts.2, aes(x = Simulations, y = Frequencies)) +
  geom_col(width = 0.6)

freq.box <- ggplot(counts.2, aes(x= Frequencies)) + geom_boxplot()
fbb <- freq.bars / freq.box
ggsave("sims_bb.png", plot = fbb, dpi = 300, scale = 2)

pf <- (pl[[1]] + pl[[2]]) / (pl[[3]] + pl[[4]]) / (pl[[5]] + pl[[6]])
ggsave("sims.png", plot = pf, dpi = 300, scale = 2)


##########################################################
##########################################################
####  DIFFERENT SIMULATIONS IN CREATE OBJECT COMMAND  ####
##########################################################
##########################################################

param <- ExtractParameters(obj1, 'All')
pv <- vector('list', length(param))
for (i in 1:length(param)){ 
  pv[[i]] <- unlist(param[[i]])
}
names(pv) <- names(param)
pv

