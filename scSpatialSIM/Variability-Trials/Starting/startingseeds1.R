library(scSpatialSIM)
library(ggplot2)
library(ggplotify)
library(patchwork)


s <- c(100, 200, 300, 400, 500, 600)
objs <- vector('list', 6)
pl <- vector('list', 6)
for (i in 1:length(s)) {
  set.seed(s[i])
  custom_window <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  # create obj
  objs[[i]] <- CreateSimulationObject(sims = 1, cell_types = 1, window = custom_window)
  
  # generate spatial pattern
  objs[[i]] <- GenerateSpatialPattern(objs[[i]], 25)
  
  # generate tissue
  objs[[i]] <- GenerateTissue(objs[[i]], density_heatmap = T, step_size = 0.1, cores = 1)

  objs[[i]] <- GenerateCellPositivity(objs[[i]], k = 1,
                                      sdmin = 0.5, sdmax = 2,
                                      density_heatmap = T, step_size = 0.1, cores = 1, probs = c(0.0, 1),
                                      shift = 1)  
  pl[[i]] <- as.ggplot(PlotSimulation(objs[[i]], which = 1, what = "whole core"))
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
  ca[[i]] <- CreateSpatialList(objs[[i]])[[1]][3:4]
  cfreq[i] <- sum(ca[[i]][2] == 1 & ca[[i]][1] == 'Tissue 2' )/sum(ca[[i]][1] == 'Tissue 2')
  cpos[i] <- sum(ca[[i]][2] == 1 & ca[[i]][1] == 'Tissue 2' )
}


f1 <- unlist(cfreq)
p1 <- unlist(cpos)
v.f1 <- var(f1)
v.p1 <- var(p1)

f1 <- c(f1, v.f1)
p1 <- c(p1, v.p1)
rs <- list(paste0('Simulation ', 1:length(objs)))[[1]]
rs <- c(rs, 'VARIANCE')
s <- c(s, NA)

counts <- data.frame(Positive = p1, Simulations= rs, Frequencies = f1, SeedValue = s)
counts

counts.2 <- counts[-(nrow(counts)), ]
freq.bars <- ggplot(counts.2, aes(x = Simulations, y = Frequencies)) +
  geom_col(width = 0.6)

freq.box <- ggplot(counts.2, aes(x= Frequencies)) + geom_boxplot()
fbb <- freq.bars / freq.box
ggsave("seeds_uncentered_bb.png", plot = fbb, dpi = 300, scale = 2)
freq.bars.s1 <- freq.bars
freq.box.s1 <- freq.box 

pf <- (pl[[1]] + pl[[2]]) / (pl[[3]] + pl[[4]]) / (pl[[5]] + pl[[6]])
ggsave("seeds_uncentered.png", plot = pf, dpi = 300, scale = 2)


ps <- s[-length(s)]
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
ps


