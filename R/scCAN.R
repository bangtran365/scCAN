#' @import scDHA
#' @importFrom FNN knn
#' @importFrom purrr quietly
#' @importFrom stats dnorm
#' @title scCAN
#' @description This is the main function to perform sc-RNA seq data clustering clustering. scCAN is fully unsupervised scRNA-seq
#' clustering framework that uses deep neural network and network fusion-based clustering algorithm.
#' First, scCAN applies a non-negative autoencoder to filter scRNA-seq data.
#' Second, the filtered data is passed to stacked Bayesian autoencoder to get multiple low-dimensional representations of input data.
#' Subsequently, scCAN converts these compressed data into networks and unify those networks to a single graph.
#' Then, scCAN uses a spectral clustering algorithm to obtain final clusters assignment.
#' @param data Gene expression matrix, with rows represent samples and columns represent genes.
#' @param sparse Boolen variable indicating whether data is a sparse matrix. The input must be a non negative sparse matrix.
#' @param n.neighbors Number of neighboring cells that are used to caculate the edge's weight. The number of neighbors are set \code{n.neighbors = 30} by default.
#' @param alpha A hyper-parameter that is used to calculate the network kernel. The value is set to \code{alpha = 0.5} by default.
#' @param n.iters A hyper-parameter to set the number of network fusion iterations. It is set to \code{n.iters = 10} by default.
#' @param ncores Number of processor cores to use.
#' @param r.seed A parameter to set a seed for reproducibility. This values is set to \code{r.seed = 1} by default.
#' @return List with the following keys:
#' \itemize{
#' \item cluster - A numeric vector containing cluster assignment for each sample. If \code{do.clus = False}, this values is always \code{NULL}.
#' \item latent - A matrix representing compressed data from the input data, with rows represent samples and columns represent latent variables.
#' }
#' @export

scCAN <- function(data, sparse = F, n.neighbors = 30, alpha = 0.5, n.iters = 10, ncores = 10L, r.seed = 1){
  set.seed(r.seed)
  if(sparse==F){
    if(max(data)>100) data <- log2(data + 1)
  }else{
    if(max(data@x)>100) data@x<- log2(data@x + 1)
  }
  result <- purrr::quietly(scDHA)(data,sparse = sparse, ncores = ncores, seed = r.seed)$result
  if(nrow(data)>5000){
    res <- cluster.big(result)
  }else{
    res <- cluster.small(result)
  }
  res
}

cluster.big <- function(data, samp.size = 5000, n.neighbors = 30, alpha = 0.5, n.iters = 10, r.seed = 1){
  set.seed(r.seed)
  K = n.neighbors
  alpha = alpha
  T = n.iters
  n_samples <- nrow(data$all.latent[[1]])
  groups <- rep(0, n_samples)
  message(n_samples)
  ind <- sample.int(n_samples, size = samp.size)

  all.sim <- list()
  for(i in 1 : length(data$all.latent)){
    dat <- data$all.latent[[i]][ind,]
    tmp <- dist2(dat,dat)^(1/2)
    all.sim[[i]] <- tmp
  }

  all.aff <- list()
  for(i in 1 : length(all.sim)){
    dat <- all.sim[[i]]
    tmp <- affinityMatrix(dat)
    all.aff[[i]] <- tmp
  }
  W = SNF(all.aff, K, T)
  res = estimateNumberOfClustersGivenGraph(W, NUMC=2:15)
  k1 <- res$`Eigen-gap best`
  k2 =res$`Eigen-gap 2nd best`
  k <- min(k1,k2)
  message(paste0("The optimal number of cluster is: ", k))

  cluster1 = spectralClustering(W,k)
  train <- data$latent[ind,]
  test  <- data$latent[-ind,]
  cluster2 <- FNN::knn(train, test, cluster1, k = 10, prob = FALSE)

  groups[ind] <-cluster1
  groups[-ind]<-cluster2

  list(cluster = groups,
       k = k,
       k1 = k1,
       k2 = k2)
}

cluster.small <- function(data, n.neighbors = 30, alpha = 0.5, n.iters = 10, r.seed = 1){
  set.seed(r.seed)
  K = n.neighbors
  alpha = alpha
  T = n.iters

  all.sim <- list()
  for(i in 1 : length(data$all.latent)){
    dat <- data$all.latent[[i]]
    tmp <- dist2(dat,dat)^(1/2)
    all.sim[[i]] <- tmp
  }

  all.aff <- list()
  for(i in 1 : length(all.sim)){
    dat <- all.sim[[i]]
    tmp <- affinityMatrix(dat)
    all.aff[[i]] <- tmp
  }

  W = SNF(all.aff, K, T)
  res = estimateNumberOfClustersGivenGraph(W, NUMC=2:15)
  k1 <- res$`Eigen-gap best`
  k2 =res$`Eigen-gap 2nd best`
  k <- min(k1,k2)
  message(paste0("The optimal number of cluster is: ", k))
  groups = spectralClustering(W,k)

  list(cluster = groups,
       k = k,
       k1 = k1,
       k2 = k2)
}

#' @title SCE
#'
#' @description SCE dataset includes scRNA-seq data and cell type information.
"SCE"
