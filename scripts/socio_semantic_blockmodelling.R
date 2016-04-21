require(blockmodeling)
require(sna)
require(igraph)
require(fpc)

source("semantic_similarity.R")
source("create_user_model.R")



getUserSimilarityMatrix <- function(userModels, threadSimMatrix, c1, c2, c3, inactives=NULL) {
  
  simMat <- matrix(1, nrow=length(userModels), ncol=length(userModels))
  
  userIds <- c()
  #maxReach <- getMaxInOutReach(userModels)
  
  usersToCheck <- which(!(1:length(userModels) %in% inactives)) # only similarities of active users have to be computed. Similarity is 1 for all others.
  print(paste("getSimilarityMatrix - ", length(userModels), "users", length(usersToCheck), "active."))
  sapply(1:(length(usersToCheck) - 1), function(a) {
    i <- usersToCheck[a]
    
    sapply((a + 1):length(usersToCheck), function(b) {
      
      j <- usersToCheck[b]
      
      forumVec1 <- append(userModels[[i]]$forums_hs, userModels[[i]]$forums_hg)
      forumVec2 <- append(userModels[[j]]$forums_hs, userModels[[j]]$forums_hg)
      
      forumSim <- 1 #getSimilaritySubforumVector(forumVec1, forumVec2)   !!! Only to speed up the calculation.
      #wordSim <- getSimilarityWordVector(wordVec1, wordVec2)
      wordSim <- (getSimilarityByMatrix(userModels[[i]]$threads_hs, userModels[[j]]$threads_hs, threadSimMatrix)
                  + getSimilarityByMatrix(userModels[[i]]$threads_hg, userModels[[j]]$threads_hg, threadSimMatrix)) / 2
      ioSim <- 1 #getSimilarityIO(userModels[[i]]$inreach, userModels[[j]]$inreach, !!! Only to speed up the calculation.
      # userModels[[i]]$outreach, userModels[[j]]$outreach)
      #print(paste("i=", i, "j=", j, " forumSim=", forumSim, " wordSim=", wordSim, " ioSim=", ioSim))
      
      simMat[j,i] <<- (c1 * forumSim + c2 * wordSim + c3 * ioSim) / (c1 + c2 + c3)
      #simMat[j,i] <<- (c1 * forumSim + c3 * ioSim) / (c1 + c3)
    })
  })
  userIds <- sapply(userModels, function(model) model$user)
  
  rownames(simMat) <- userIds
  colnames(simMat) <- userIds
  simMat
}

getSemanticBlockmodel <- function(network, postData, threadSimMatrix, minTime, maxTime, numClusters=NULL, userSimMatrix=NULL) {
  
  if (is.null(E(network)$weight)) {
    
    E(network)$weight <- 1
  }
  
  if (is.null(userSimMatrix)) {
    
    print("create user models")
    userModels <- lapply(V(network)$label, getUserModelHSHG, postData, network, minTime, maxTime)
    save(userModels, file="user_models.Rdata")
    
    print("compute similarity matrix")
    userSimMatrix <- getUserSimilarityMatrix(userModels, threadSimMatrix, 1, 1, 1)
    save(userSimMatrix, file="similarities.Rdata")  
  }
  
  print("fit blockmodel")
  d <- as.dist(1 - userSimMatrix)
  h <- hclust(d, method = "ward.D")
  
  if (is.null(numClusters)) {
    
    numClusters <- nselectboot(d, clustermethod=disthclustCBI, method="ward.D")$kopt
    print(paste("num. clusters", numClusters))
  }
  
  part <- cutree(h, numClusters)
  #   partMapped <- sapply(V(network)$label, function(nodeLabel) {
  #     
  #     simRow <- which(rownames(sim) == nodeLabel)
  #     
  #     ifelse((length(simRow) > 0), part[simRow], (max(part) + 1))
  #   })
  
  # isolated nodes should not contribute to the denstity
  nodesConnected <- length(which(igraph::degree(network) > 0))
  dens <- length(E(network)) / (nodesConnected * (nodesConnected - 1))
  crit.fun(get.adjacency(network, sparse=FALSE), part, approach="bin", blocks=c("null","reg"), blockWeights=c(null=1,reg=dens/(1-dens)), norm=TRUE)
}


start <- function() {
  
  # Testrun
  network <- read.graph("testgraph.gml", "gml")
  postData <- read.csv2("../../network_extraction/post_tables/finance_001_featurised/finance_001_classified.csv")
  threadSimMatrix <- readMM("../finance_001/thread_sem_matrix_wordnet.txt") #read.csv("../finance_001/thread_keywords_finance_001/calais_annotations/calais_tags_finance.csv")
  dimNames <- read.table("../finance_001/thread_sem_sim_wordnet_dimnames.txt", as.is = "x")$x
  
  rownames(threadSimMatrix) <- dimNames
  colnames(threadSimMatrix) <- dimNames
  
  bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, 10)
  
  save(bm, file="blockmodel.Rdata")
  
  bm
  
  # TODO: 
  ## More classified datasets
  ## Try thread similartiy based on LSA and compare with wordnet approach
  ## Determine time window size
  ## Information diffusion by IM_slice1 X IM_slice2 X ... (implemented)
  ## Combine semantic and structural blockmodels (implemented)
}