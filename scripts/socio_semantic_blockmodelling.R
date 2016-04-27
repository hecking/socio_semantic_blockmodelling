require(blockmodeling)
require(sna)
require(igraph)
require(fpc)
require(psych)

getSocioSemanticBlockmodel <- function(network, nonNetworkSimilarity, sigmaNetSim, sigmaNonNetSim, numClusters=NULL) {
  
  if (is.null(E(network)$weight)) {
    
    E(network)$weight <- 1
  }
  
  userRegeSim <- REGE.for(get.adjacency(network, sparse = FALSE))$E
  userSim <- (sigmaNetSim * userRegeSim + sigmaNonNetSim * nonNetworkSimilarity) / (sigmaNetSim + sigmaNonNetsim)
  
  print("fit blockmodel")
  d <- as.dist(1 - userSim)
  h <- hclust(d, method = "ward.D")
  
  if (is.null(numClusters)) {
    
    numClusters <- nselectboot(d, clustermethod=disthclustCBI, method="ward.D")$kopt
    print(paste("num. clusters", numClusters))
  }
  
  part <- cutree(h, numClusters)
  
  nodesConnected <- length(which(igraph::degree(network) > 0))
  dens <- length(E(network)) / (nodesConnected * (nodesConnected - 1))
  crit.fun(get.adjacency(network, sparse=FALSE), part, approach="bin", blocks=c("null","reg"), blockWeights=c(null=1,reg=dens/(1-dens)), norm=TRUE)
}

getSimilarityCorrelations <- function(structuralSimilarity, regularSimilarity, semanticSimilarity) {
  
  data <- data.frame(
    "structural"=as.vector(structuralSimilarity[lower.tri(structuralSimilarity)]),
    "regular"=as.vector(regularSimilarity[lower.tri(regularSimilarity)]),
    "semantic"=as.vector(semanticSimilarity[lower.tri(semanticSimilarity)])
  )
  
  corr.test(data, method="spearman")
}

start <- function() {
  
  network <- read.graph("example_net.gml", "gml")
  toDelete <- which(igraph::degree(network) == 0)
  network <- delete.vertices(network, toDelete)
  userSemanticSimilarity <- read.csv("example_sem_sim.csv")
  userRegularSimililarity <- REGE.for(get.adjacency(network, sparse=FALSE))$E
  userStructuralSimilarity1 <- 1 - sedist(get.adjacency(network, sparse=FALSE), method="euclidean")  
  
  # similarity correlations
  simCor <- getSimilarityCorrelations(userStrSim, userRegSim, userSemSim)
  print("similarity correlations")
  print("Spearman correlations")
  print(simCor$r)
  print("p values")
  print(simCor$p)
  
  # --------------------------------------------------
  # semantic similarity = 2, regular similarity = 1
  # --------------------------------------------------
  
  bm <- getSocioSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=userRegSim)
  
  save(bm, file="example.Rdata")
  
  bm
}