spatial_data = CreateSimulationObject(sims = 1, cell_types = 1) %>%
  #produce the point pattern
  GenerateSpatialPattern() %>%
  #make tissues
  GenerateTissue(density_heatmap = FALSE, step_size = 0.1, cores = 1) %>%
  #create positive and negative cells
  GenerateCellPositivity(k = 4, sdmin = 3, sdmax = 5,
                         density_heatmap = FALSE, step_size = 1, cores = 1, probs = c(0.0, 0.1), shift = 0) %>%
  #convert to a list of spatial data frames
  CreateSpatialList(single_df = FALSE)
spat_data_distribution = GenerateDistributions(spatial_data)



spat_data_distribution


CreateSpatialList(sim_obj_1)
