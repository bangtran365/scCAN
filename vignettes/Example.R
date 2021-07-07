## ---- eval=FALSE--------------------------------------------------------------
#  library(scCAN)
#  #Load example data (SCE dataset)
#  data("SCE")
#  
#  #Get data matrix and label
#  data <- t(SCE$data); label <- as.character(SCE$cell_type1)
#  

## ---- eval=FALSE--------------------------------------------------------------
#  #Generate clustering result, the input matrix has rows as samples and columns as genes
#  result <- scCAN(data, r.seed = 1)
#  
#  #The clustering result can be found here
#  cluster <- result$cluster
#  
#  #Calculate adjusted Rand Index using mclust package
#  ari <- round(scCAN::adjustedRandIndex(cluster,label), 2)
#  print(paste0("ARI = ", ari))

