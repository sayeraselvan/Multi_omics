library(kuenm)

variables <- "tif"  

kuenm_rpca(
  variables = variables, 
  in.format = "GTiff", 
  var.scale = TRUE, 
  write.result = TRUE, 
  out.format = "ascii", 
  out.dir = "PCA_results", 
  project = FALSE
)

