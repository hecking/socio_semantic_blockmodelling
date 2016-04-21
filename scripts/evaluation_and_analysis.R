require(blockmodeling)
require(sna)
require(igraph)
require(psych)

source("create_user_model.R")

getSimilarityCorrelations <- function(structuralSimilarity, regularSimilarity, semanticSimilarity) {
  
  data <- data.frame(
    "structural"=as.vector(structuralSimilarity[lower.tri(structuralSimilarity)]),
    "regular"=as.vector(regularSimilarity[lower.tri(regularSimilarity)]),
    "semantic"=as.vector(semanticSimilarity[lower.tri(semanticSimilarity)])
  )
  
  corr.test(data, method="spearman")
}

evaluateBlockmodel <- function(userModels, threadKeywords, blockmodel, userDistances, filePrefix) {
  
  print("evaluate blockmodel")
  # for each cluster
  write.csv(blockmodel$clu, paste(filePrefix, "membership.csv", sep="_"))  
  res <- t(sapply(unique(blockmodel$clu), function(cid) {
    
    members <- which(blockmodel$clu == cid)
    sublist <- userModels[members]
    
    outreach <- 0
    inreach <- 0
    wordsHG <- c()
    forumsHG <- c()
    wordsHS <- c()
    forumsHS <- c()
    
    sapply(sublist, function(userModel) {
      
      outreach <<- outreach + userModel$outreach
      inreach <<- inreach + userModel$inreach
      forumsHG <<- append(forumsHG, userModel$forums_hg)
      forumsHS <<- append(forumsHS, userModel$forums_hs)
      
      keywordsHS <- threadKeywords$keyword[which(threadKeywords$thread_id %in% userModel$threads_hs)]
      keywordsHG <- threadKeywords$keyword[which(threadKeywords$thread_id %in% userModel$threads_hg)]
      
      if (length(keywordsHS) > 0) {
        
        wordsHS <<- append(wordsHS, keywordsHS)
      }
      if (length(keywordsHG) > 0) {
        
        wordsHG <<- append(wordsHG, keywordsHG)
      }
    })
    
    ## get mean outreach
    ## get mean inreach
    ## get 5 most frequent words ranked
    wordFreqHS <- sort(table(wordsHS), decreasing = TRUE)
    wordFreqHG <- sort(table(wordsHG), decreasing = TRUE)
    
    ## get 5 most frequent forums ranked
    forumFreqHS <- sort(table(forumsHS), decreasing = TRUE)
    forumFreqHG <- sort(table(forumsHG), decreasing = TRUE)
    
    
    c(cid, length(members), mean(outreach), mean(inreach), 
      ifelse(length(forumFreqHS) > 0, names(forumFreqHS[1]), "none"),
      ifelse(length(forumFreqHS) > 1, names(forumFreqHS[2]), "none"),
      ifelse(length(forumFreqHS) > 2, names(forumFreqHS[3]), "none"),
      ifelse(length(forumFreqHG) > 0, names(forumFreqHG[1]), "none"),
      ifelse(length(forumFreqHG) > 1, names(forumFreqHG[2]), "none"),
      ifelse(length(forumFreqHG) > 2, names(forumFreqHG[3]), "none"),
      ifelse(length(wordFreqHS) > 0, names(wordFreqHS[1]), "none"),
      ifelse(length(wordFreqHS) > 1, names(wordFreqHS[2]), "none"),
      ifelse(length(wordFreqHS) > 2, names(wordFreqHS[3]), "none"),
      ifelse(length(wordFreqHS) > 3, names(wordFreqHS[4]), "none"),
      ifelse(length(wordFreqHS) > 4, names(wordFreqHS[5]), "none"),
      ifelse(length(wordFreqHG) > 0, names(wordFreqHG[1]), "none"),
      ifelse(length(wordFreqHG) > 1, names(wordFreqHG[2]), "none"),
      ifelse(length(wordFreqHG) > 2, names(wordFreqHG[3]), "none"),
      ifelse(length(wordFreqHG) > 3, names(wordFreqHG[4]), "none"),
      ifelse(length(wordFreqHG) > 4, names(wordFreqHG[5]), "none")
    )
  }))
  
  
  # output
  ## cluster#num_members#mean_outreach#mean_inreach1stforum_freqency#2ndforum_frequency#...#1stword_frequency#2ndword_frequency#...
  colnames(res) <- c("cluster", "num_members", "mean_outreach", "mean_inreach", 
                     "hs_forum1_(freq)", "hs_forum2_(freq)", "hs_forum3_(freq)",
                     "hg_forum1_(freq)", "hg_forum2_(freq)", "hg_forum3_(freq)",
                     "hs_word1_(freq)", "hs_word2_(freq)", "hs_word3_(freq)","hs_word4_(freq)", "hs_word5_(freq)",
                     "hg_word1_(freq)", "hg_word2_(freq)", "hg_word3_(freq)","hg_word4_(freq)", "hg_word5_(freq)")
  write.csv(res, paste(filePrefix, "clusters.csv", sep="_"))
  
  print("cluster content analysis done!")
  
  ## Image matrix
  write.csv(blockmodel$IM, paste(filePrefix, "image_matrix.csv", sep="_"))
  
  ## Block error matrix
  write.csv(blockmodel$E, paste(filePrefix, "block_error_matrix.csv", sep="_"))
  
  ## Overall error
  print(blockmodel$err)
  
  ## Role interactions
  am <- matrix(FALSE, nrow=nrow(blockmodel$IM), ncol=ncol(blockmodel$IM))
  am[which(blockmodel$IM != "null")] <- TRUE
  g <- graph.adjacency(am)
  
  png(paste(filePrefix, "role_interactions.png", sep="_"))
  plot(g, vertex.size=(as.numeric(res[sort(res[,1], index.return=TRUE)$ix,2]) / 5), edge.width=(0.01 / blockmodel$E[am]))
  dev.off()
  
  print("blockmodel analysis done!")
  
  cstats <- cluster.stats(userDistances, clustering = blockmodel$clu)
  
  cstatsTable <- data.frame(cases=cstats$n, num_clusters=cstats$cluster.number, avg_inner_distance=cstats$average.within, avg_outer_distance=cstats$average.between, dunn_index=cstats$dunn, bm_error=blockmodel$err, bm_dunn_ratio=(blockmodel$err / cstats$dunn))
  
  write.csv(cstatsTable, paste(filePrefix, "cluster_statistics.csv", sep="_"))
  
  print("cluster statistics done!")
  
  cstatsTable
}

evaluateDynamicBlockmodel <- function(userModelsList, threadKeywords, blockmodels, userDistances, networks, filePrefix) {
  
  membershipVec <- blockmodels[[1]]$clu # Clusters are the same for every time slice.
  clusters <- lapply(unique(membershipVec), function(cl) {
    
    which(membershipVec == cl)
  })
  numClusters <- length(unique(membershipVec))
  evaluateBlockmodel(userModelsList[[1]], threadKeywords, blockmodels[[1]], userDistances, paste(filePrefix, "slice_1", sep="_"))
  am <- matrix(FALSE, nrow=nrow(blockmodels[[1]]$IM), ncol=ncol(blockmodels[[1]]$IM))
  am[which(blockmodels[[1]]$IM != "null")] <- TRUE
  #diffusionGraph <- graph.adjacency(am)
  diffusionGraph <- graph.full(numClusters, loops = TRUE, directed = TRUE) #graph.adjacency(blockmodels[[1]]$E, weighted = TRUE)
  #E(diffusionGraph)$time1 <- 1
  #E(diffusionGraph)$weight <- 0.01 / blockmodels[[1]]$E[am]
  E(diffusionGraph)$weight1 <- sapply(E(diffusionGraph), function(e) {ifelse(am[e], 0.01 / blockmodels[[1]]$E[e], 100 * blockmodels[[1]]$E[e])})
  
  activeUsers <- which(igraph::degree(networks[[1]]) > 0)
  
  activity <- sapply(clusters, function(cluster) {
    
    length(intersect(cluster, activeUsers)) / length(cluster)
  })
  
  if (length(blockmodels) > 1) {
    
    activity <- cbind(activity, sapply(2:length(blockmodels), function(i) {
      
      evaluateBlockmodel(userModelsList[[i]], threadKeywords, blockmodels[[i]], userDistances, paste(filePrefix, "slice", i, sep="_"))
      
      am <- matrix(FALSE, nrow=nrow(blockmodels[[i]]$IM), ncol=ncol(blockmodels[[i]]$IM))
      am[which(blockmodels[[i]]$IM != "null")] <- TRUE
      #diffusionGraphNew <- graph.adjacency(am)
      diffusionGraphNew <- graph.full(numClusters, loops = TRUE, directed = TRUE) #graph.adjacency(blockmodels[[i]]$E, weighted = TRUE)
      #E(diffusionGraphNew)$time <- i
      #E(diffusionGraphNew)$weight <- 0.01 / blockmodels[[i]]$E[am]
      eWeights <- sapply(E(diffusionGraphNew), function(e) {ifelse(am[e], 0.01 / blockmodels[[i]]$E[e], 100 * blockmodels[[i]]$E[e])})
      diffusionGraphNew <- set.edge.attribute(diffusionGraphNew, paste("weight", i, sep=""), value=eWeights)
      diffusionGraph <<- diffusionGraph %u% diffusionGraphNew
      
      activeUsers <- which(igraph::degree(networks[[i]]) > 0)
      
      sapply(clusters, function(cluster) {
        
        length(intersect(cluster, activeUsers)) / length(cluster)
      })
    }))
  }
  
  write.graph(diffusionGraph, paste(filePrefix, "diffusion_graph.gml", sep="_"), "gml")
  
  colnames(activity) <- paste("slice", 1:length(blockmodels), sep="_")
  rownames(activity) <- as.character(unique(clusters))
  
  write.csv(activity, paste(filePrefix, "cluster_activity.csv", sep="_"))
  
  activity
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

### OLD DATA ###

# evaluate <- function() {
#   
#   postData <- read.csv2("../../network_extraction/post_tables/finance_001_featurised/finance_001_classified.csv")
#   keywords <- read.csv("../finance_001/thread_keywords_finance_001/calais_annotations/calais_tags_finance.csv", as.is="keyword")
#   threadSimMatrix <- readMM("../finance_001/thread_sem_matrix_wordnet.txt")
#   dimNames <- read.table("../finance_001/thread_sem_sim_wordnet_dimnames.txt", as.is = "x")$x
#   rownames(threadSimMatrix) <- dimNames
#   colnames(threadSimMatrix) <- dimNames
#   
# #   threadSimMatrix <- readMM("../finance_001/thread_sem_matrix_lsa.txt")
# #   dimNames <- read.table("../finance_001/thread_sem_sim_lsa_dimnames.txt", as.is = "x")$x
# #   rownames(threadSimMatrix) <- dimNames
# #   colnames(threadSimMatrix) <- dimNames
#   
#   graphDir <- "../../network_extraction/classified_datasets/fincance_001/timeslices_35_86_w3/"
#   graphFiles <- list.files(graphDir, full.names = TRUE)
#   graphs <- lapply(graphFiles, function(file) {
#     
#     read.graph(file, "gml")
#   })
#   
#   network <- collapseGraphs(graphs)
#   toDelete <- which(igraph::degree(network) == 0)
#   network <- delete.vertices(network, toDelete)
#   
#   outputFilePrefix <- "finance_001"
#   
#   ####################################################
#   #
#   # Static blockmodels
#   #
#   ####################################################
#   
# #   minTime <- 35
# #   maxTime <- 86
# #   print("user models static")
# #   userModels <- lapply(V(network)$label, getUserModelHSHG, postData, network, minTime, maxTime)
# #   save(userModels, file="user_model_static.Rdata")
# #   print("user similarity static")
# #   userSemSim <- getUserSimilarityMatrix(userModels, threadSimMatrix, 1, 1, 1)
# #   save(userSemSim, file="user_sim_static.Rdata")
# #   userRegSim <- REGE.for(get.adjacency(network, sparse=FALSE))$E
# #   userStrSim <- 1 - sedist(get.adjacency(network, sparse=FALSE), method="euclidean")
# #   
# #   # similarity correlations
# #   simCor <- getSimilarityCorrelations(userStrSim, userRegSim, userSemSim)
# #   write.csv(simCor$R, paste("results/static_similarity_correlations_", outputFilePrefix, ".csv", sep=""))
# #   write.csv(simCor$P, paste("results/static_similarity_correlations_p_", outputFilePrefix, ".csv", sep=""))
# #   
# #   
# #   # --------------------------------------------------
# #   # semantic similarity = 0, regular similarity = 1
# #   # --------------------------------------------------
# #   
# #   bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=userRegSim)
# #   distances <- as.dist(1 - userRegSim)
# #   
# #   bmstats_static_s0_r1 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results/static_s0_r1", outputFilePrefix, sep="_"))
# #   
# #   
# #   # --------------------------------------------------
# #   # semantic similarity = 1, regular similarity = 0
# #   # --------------------------------------------------
# #   
# #   bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=userSemSim)
# #   distances <- as.dist(1 - userSemSim)
# #   
# #   bmstats_static_s1_r0 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results/static_s1_r0", outputFilePrefix, sep="_"))
# #   
# #   
# #   # --------------------------------------------------
# #   # semantic similarity = 1, regular similarity = 1
# #   # --------------------------------------------------
# #   
# #   sim <- (userSemSim + userRegSim) / 2
# #   bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=sim)
# #   distances <- as.dist(1 - sim)
# #   
# #   bmstats_static_s1_r1 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results/static_s1_r1", outputFilePrefix, sep="_"))
# #   
# #   
# #   # --------------------------------------------------
# #   # semantic similarity = 2, regular similarity = 1
# #   # --------------------------------------------------
# # 
# #   sim <- (2 * userSemSim + userRegSim) / 3
# #   bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=sim)
# #   distances <- as.dist(1 - sim)
# #   
# #   bmstats_static_s2_r1 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results/static_s2_r1", outputFilePrefix, sep="_"))
# #   
# #   
# #   # --------------------------------------------------
# #   # semantic similarity = 1, regular similarity = 2
# #   # --------------------------------------------------
# #   
# #   sim <- (userSemSim + 2 * userRegSim) / 3
# #   bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=sim)
# #   distances <- as.dist(1 - sim)
# #   
# #   bmstats_static_s1_r2 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results/static_s1_r2", outputFilePrefix, sep="_"))
# #   
# #   
# #   #-------
# #   
# #   bmStatsTable <- rbind(bmstats_static_s0_r1, bmstats_static_s1_r2, bmstats_static_s1_r1, bmstats_static_s2_r1, bmstats_static_s1_r0)
# #   
# #   write.csv(paste(bmStatsTable, "results/bmstats_static", outputFilePrefix, ".csv", sep=""))
# #   #-------
#   
#   
#   ####################################################
#   #
#   # Dynamic blockmodels
#   #
#   ####################################################
#   
#    timeSliceSize <- floor(length(graphs) / 3)
#    
#   network1 <- collapseGraphs(graphs[1:timeSliceSize])
#   network1 <- delete.vertices(network1, toDelete)
#   network2 <- collapseGraphs(graphs[(timeSliceSize + 1):(2 * timeSliceSize)])
#   network2 <- delete.vertices(network2, toDelete)
#   network3 <- collapseGraphs(graphs[(2 * timeSliceSize + 1):length(graphs)])
#   network3 <- delete.vertices(network3, toDelete)
#   
# #   print("user models dynamic1")
# #   minTime <- 35
# #   maxTime <- 49
# #   inactives <- which(igraph::degree(network1) == 0)
# #   userModels1 <- lapply(V(network)$label, getUserModelHSHG, postData, network1, minTime, maxTime)
# #   save(userModels1, file="user_model_dynamic1.Rdata")
#   print("user similarity dynamic1")
# #   userSemSimilarity1 <- getUserSimilarityMatrix(userModels1, threadSimMatrix, 1, 1, 1, inactives)
# #   save(userSemSimilarity1, file="user_sim_dynamic1.Rdata")
#   userRegSimilarity1 <- REGE.for(get.adjacency(network1, sparse=FALSE))$E
#   userStrSimilarity1 <- 1 - sedist(get.adjacency(network1, sparse=FALSE), method="euclidean")
#   
# #   print("user models dynamic2")
# #   minTime <- 50
# #   maxTime <- 64
# #   inactives <- which(igraph::degree(network2) == 0)
# #   userModels2 <- lapply(V(network)$label, getUserModelHSHG, postData, network2, minTime, maxTime)
# #   save(userModels2, file="user_model_dynamic2.Rdata")
#   print("user similarity dynamic2")
# #   userSemSimilarity2 <- getUserSimilarityMatrix(userModels2, threadSimMatrix, 1, 1, 1, inactives)
# #   save(userSemSimilarity2, file="user_sim_dynamic2.Rdata")
#   userRegSimilarity2 <- REGE.for(get.adjacency(network2, sparse=FALSE))$E
#   userStrSimilarity2 <- 1 - sedist(get.adjacency(network2, sparse=FALSE), method="euclidean")
#    
# #  print("user models dynamic3")
# #  minTime <- 65
# #  maxTime <- 86
# #  inactives <- which(igraph::degree(network3) == 0)
# #   userModels3 <- lapply(V(network)$label, getUserModelHSHG, postData, network3, minTime, maxTime)
# #   save(userModels3, file="user_model_dynamic3.Rdata")
#   print("user similarity dynamic3")
# #   userSemSimilarity3 <- getUserSimilarityMatrix(userModels3, threadSimMatrix, 1, 1, 1, inactives)
# #   save(userSemSimilarity3, file="user_sim_dynamic3.Rdata")
#   userRegSimilarity3 <- REGE.for(get.adjacency(network3, sparse=FALSE))$E
#   userStrSimilarity3 <- 1 - sedist(get.adjacency(network3, sparse=FALSE), method="euclidean")
# 
# 
#   userSemSim <- userSemSimilarity1 + userSemSimilarity2 + userSemSimilarity3
#   userRegSim <- userRegSimilarity1 + userRegSimilarity2 + userRegSimilarity3
#   userStrSim <- userStrSimilarity1 + userStrSimilarity2 + userStrSimilarity3
#   
#   
#   # similarity correlations
# #   simCor <- getSimilarityCorrelations(userStrSim, userRegSim, userSemSim)
# #   write.csv(simCor$R, paste("results/dynamic_similarity_correlations_", outputFilePrefix, ".csv", sep=""))
# #   write.csv(simCor$P, paste("results/dynamic_similarity_correlations_p_", outputFilePrefix, ".csv", sep=""))
#   
#   
#   # --------------------------------------------------
#   # semantic similarity = 0, regular similarity = 1
#   # --------------------------------------------------
#   
#   distances <- as.dist(1 - userRegSim)
#   nc <- nselectboot(distances, clustermethod=disthclustCBI, method="ward.D")$kopt
#   print(paste("num clusters:", nc))
#   if (nc < 2) {
#     
#     nc <- 2
#   }
#   bm1 <- getSemanticBlockmodel(network1, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=userRegSim)
#   bm2 <- getSemanticBlockmodel(network2, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=userRegSim)
#   bm3 <- getSemanticBlockmodel(network3, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=userRegSim)
#   
#   bmstats_dynamic_s0_r1 <- evaluateDynamicBlockmodel(list(userModels1, userModels2, userModels3), keywords, list(bm1, bm2, bm3), distances, 
#                                                      list(network1, network2, network3), paste("results/dynamic_s0_r1", outputFilePrefix, sep="_"))
#   
#   
#   # --------------------------------------------------
#   # semantic similarity = 1, regular similarity = 0
#   # --------------------------------------------------
#   
#   distances <- as.dist(1 - userSemSim)
#   nc <- nselectboot(distances, clustermethod=disthclustCBI, method="ward.D")$kopt
#   
#   if (nc < 2) {
#     
#     nc <- 2
#   }
# 
#   bm1 <- getSemanticBlockmodel(network1, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=userSemSim)
#   bm2 <- getSemanticBlockmodel(network2, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=userSemSim)
#   bm3 <- getSemanticBlockmodel(network3, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=userSemSim)
#   
#   bmstats_dynamic_s1_r0 <- evaluateDynamicBlockmodel(list(userModels1, userModels2, userModels3), keywords, list(bm1, bm2, bm3), distances, 
#                                                      list(network1, network2, network3), paste("results/dynamic_s1_r0", outputFilePrefix, sep="_"))
#   
#   
#   # --------------------------------------------------
#   # semantic similarity = 1, regular similarity = 1
#   # --------------------------------------------------
#   
#   sim <- (userSemSim + userRegSim) / 2
#   
#   distances <- as.dist(1 - sim)
#   nc <- nselectboot(distances, clustermethod=disthclustCBI, method="ward.D")$kopt
# 
#   if (nc < 2) {
#     
#     nc <- 2
#   }
# 
#   bm1 <- getSemanticBlockmodel(network1, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=sim)
#   bm2 <- getSemanticBlockmodel(network2, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=sim)
#   bm3 <- getSemanticBlockmodel(network3, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=sim)
#   
#   bmstats_dynamic_s1_r1 <- evaluateDynamicBlockmodel(list(userModels1, userModels2, userModels3), keywords, list(bm1, bm2, bm3), distances, 
#                                                      list(network1, network2, network3), paste("results/dynamic_s1_r1", outputFilePrefix, sep="_"))
#   
#   
#   # --------------------------------------------------
#   # semantic similarity = 2, regular similarity = 1
#   # --------------------------------------------------
#   
#   sim <- (2 * userSemSim + userRegSim) / 3
#   
#   distances <- as.dist(1 - sim)
#   nc <- nselectboot(distances, clustermethod=disthclustCBI, method="ward.D")$kopt
# 
#   if (nc < 2) {
#     
#     nc <- 2
#   }
# 
#   bm1 <- getSemanticBlockmodel(network1, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=sim)
#   bm2 <- getSemanticBlockmodel(network2, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=sim)
#   bm3 <- getSemanticBlockmodel(network3, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=sim)
#   
#   bmstats_dynamic_s2_r1 <- evaluateDynamicBlockmodel(list(userModels1, userModels2, userModels3), keywords, list(bm1, bm2, bm3), distances, 
#                                                      list(network1, network2, network3), paste("results/dynamic_s2_r1", outputFilePrefix, sep="_"))
#   
#   
#   # --------------------------------------------------
#   # semantic similarity = 1, regular similarity = 2
#   # --------------------------------------------------
#   
#   sim <- (userSemSim + 2 * userRegSim) / 3
#   
#   distances <- as.dist(1 - sim)
#   nc <- nselectboot(distances, clustermethod=disthclustCBI, method="ward.D")$kopt
#   
#   if (nc < 2) {
#     
#     nc <- 2
#   }
# 
#   bm1 <- getSemanticBlockmodel(network1, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=sim)
#   bm2 <- getSemanticBlockmodel(network2, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=sim)
#   bm3 <- getSemanticBlockmodel(network3, postData, threadSimMatrix, minTime, maxTime, numClusters = nc, userSimMatrix=sim)
#   
#   bmstats_dynamic_s1_r2 <- evaluateDynamicBlockmodel(list(userModels1, userModels2, userModels3), keywords, list(bm1, bm2, bm3), distances, 
#                                                      list(network1, network2, network3), paste("results/dynamic_s1_r2", outputFilePrefix, sep="_"))
#   
#   
#   #-------
#     
#   bmStatsTable <- rbind(bmstats_dynamic_s0_r1, bmstats_dynamic_s1_r2, bmstats_dynamic_s1_r1, bmstats_dynamic_s2_r1, bmstats_dynamic_s1_r0)
#   
#   write.csv(paste(bmStatsTable, "results/activities_dynamic", outputFilePrefix, ".csv", sep=""))
#   #-------
# }
evaluate <- function() {
  
  postData <- read.csv2("../../network_extraction/post_tables/finance_001_featurised/latest/classified.csv")
  keywords <- read.csv("../finance_001/thread_keywords_finance_001/calais_annotations/calais_tags_finance.csv", as.is="keyword")
  threadSimMatrix <- readMM("../finance_001/thread_sem_matrix_wordnet.txt")
  dimNames <- read.table("../finance_001/thread_sem_sim_wordnet_dimnames.txt", as.is = "x")$x
  rownames(threadSimMatrix) <- dimNames
  colnames(threadSimMatrix) <- dimNames
  threadSimMatrix[which(threadSimMatrix > 1)] <- threadSimMatrix[which(threadSimMatrix > 1)] / 1e+308
  
  #   threadSimMatrix <- readMM("../finance_001/thread_sem_matrix_lsa.txt")
  #   dimNames <- read.table("../finance_001/thread_sem_sim_lsa_dimnames.txt", as.is = "x")$x
  #   rownames(threadSimMatrix) <- dimNames
  #   colnames(threadSimMatrix) <- dimNames
  
  graphDir <- "../../network_extraction/classified_datasets/fincance_001/latest/timeslices/"
  graphFiles <- list.files(graphDir, full.names = TRUE)
  graphs <- lapply(graphFiles, function(file) {
    
    read.graph(file, "gml")
  })
  
  network <- collapseGraphs(graphs)
  toDelete <- which(igraph::degree(network) == 0)
  network <- delete.vertices(network, toDelete)
  
  outputFilePrefix <- "finance_001"
  
  ####################################################
  #
  # Static blockmodels
  #
  ####################################################
  
  minTime <- 1
  maxTime <- 70
  print("user models static")
  userModels <- lapply(V(network)$label, getUserModelHSHG, postData, network, minTime, maxTime)
  #save(userModels, file="user_model_static.Rdata")
  print("user similarity static")
  #userSemSim <- getUserSimilarityMatrix(userModels, threadSimMatrix, 0, 1, 0)
  #save(userSemSim, file="results_finance_001/user_sim_static.Rdata")
  load("results_finance_001/user_sim_static.Rdata")
  userSemSim[which(userSemSim == Inf)] <- 0
  userRegSim <- REGE.for(get.adjacency(network, sparse=FALSE))$E
  userStrSim <- 1 - sedist(get.adjacency(network, sparse=FALSE), method="euclidean")
  
  # similarity correlations
  simCor <- getSimilarityCorrelations(userStrSim, userRegSim, userSemSim)
  write.csv(simCor$r, paste("results_finance_001/static_similarity_correlations_", outputFilePrefix, ".csv", sep=""))
  write.csv(simCor$p, paste("results_finance_001/static_similarity_correlations_p_", outputFilePrefix, ".csv", sep=""))
  
  distances <- as.dist(1 - userSemSim)
  # --------------------------------------------------
  # semantic similarity = 0, regular similarity = 1
  # --------------------------------------------------
  
  bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=userRegSim)
  #distances <- as.dist(1 - userRegSim)
  
  bmstats_static_s0_r1 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results_finance_001/static_s0_r1", outputFilePrefix, sep="_"))
  
  
  # --------------------------------------------------
  # semantic similarity = 1, regular similarity = 0
  # --------------------------------------------------
  
  bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=userSemSim)
  #distances <- as.dist(1 - userSemSim)
  
  bmstats_static_s1_r0 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results_finance_001/static_s1_r0", outputFilePrefix, sep="_"))
  
  
  # --------------------------------------------------
  # semantic similarity = 1, regular similarity = 1
  # --------------------------------------------------
  
  sim <- (userSemSim + userRegSim) / 2
  bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=sim)
  #distances <- as.dist(1 - sim)
  
  bmstats_static_s1_r1 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results_finance_001/static_s1_r1", outputFilePrefix, sep="_"))
  
  
  # --------------------------------------------------
  # semantic similarity = 2, regular similarity = 1
  # --------------------------------------------------
  
  sim <- (2 * userSemSim + userRegSim) / 3
  bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=sim)
  #distances <- as.dist(1 - sim)
  
  bmstats_static_s2_r1 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results_finance_001/static_s2_r1", outputFilePrefix, sep="_"))
  
  # --------------------------------------------------
  # semantic similarity = 3, regular similarity = 1
  # --------------------------------------------------
  
  sim <- (3 * userSemSim + userRegSim) / 4
  bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=sim)
  #distances <- as.dist(1 - sim)
  
  bmstats_static_s3_r1 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results_finance_001/static_s3_r1", outputFilePrefix, sep="_"))
  
  # --------------------------------------------------
  # semantic similarity = 4, regular similarity = 1
  # --------------------------------------------------
  
  sim <- (4 * userSemSim + userRegSim) / 5
  bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=sim)
  #distances <- as.dist(1 - sim)
  
  bmstats_static_s4_r1 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results_finance_001/static_s4_r1", outputFilePrefix, sep="_"))
  # --------------------------------------------------
  # semantic similarity = 1, regular similarity = 2
  # --------------------------------------------------
  
  sim <- (userSemSim + 2 * userRegSim) / 3
  bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=sim)
  #distances <- as.dist(1 - sim)
  
  bmstats_static_s1_r2 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results_finance_001/static_s1_r2", outputFilePrefix, sep="_"))
  
  # --------------------------------------------------
  # semantic similarity = 1, regular similarity = 3
  # --------------------------------------------------
  
  sim <- (userSemSim + 3 * userRegSim) / 4
  bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=sim)
  #distances <- as.dist(1 - sim)
  
  bmstats_static_s1_r3 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results_finance_001/static_s1_r3", outputFilePrefix, sep="_"))
  
  # --------------------------------------------------
  # semantic similarity = 1, regular similarity = 4
  # --------------------------------------------------
  
  sim <- (userSemSim + 4 * userRegSim) / 5
  bm <- getSemanticBlockmodel(network, postData, threadSimMatrix, 0, 0, userSimMatrix=sim)
  #distances <- as.dist(1 - sim)
  
  bmstats_static_s1_r4 <- evaluateBlockmodel(userModels, keywords, bm, distances, paste("results_finance_001/static_s1_r4", outputFilePrefix, sep="_"))
  
  #-------
  
  bmStatsTable <- rbind(bmstats_static_s0_r1, bmstats_static_s1_r4, bmstats_static_s1_r3, bmstats_static_s1_r2, bmstats_static_s1_r1, bmstats_static_s2_r1, bmstats_static_s3_r1, bmstats_static_s4_r1, bmstats_static_s1_r0)
  
  row.names(bmStatsTable) <- c("sem_only", "s4_r1", "s3_r1", "s2_r1", "s1_r1", "s1_r2", "s1_r3", "s1_r4", "reg_only")
  write.csv(bmStatsTable, paste("results_finance_001/bmstats_static_", outputFilePrefix, ".csv", sep=""))
  #-------
  
  #   
  #   ####################################################
  #   #
  #   # Dynamic blockmodels
  #   #
  #   ####################################################
  #   
  
  network1 <- collapseGraphs(graphs[1:2])
  toDelete <- which(igraph::degree(network1) == 0)
  network1 <- delete.vertices(network1, toDelete)
  
  network2 <- collapseGraphs(graphs[3:4])
  toDelete <- which(igraph::degree(network2) == 0)
  network2 <- delete.vertices(network2, toDelete)
  
  network3 <- collapseGraphs(graphs[5:6])
  toDelete <- which(igraph::degree(network3) == 0)
  network3 <- delete.vertices(network3, toDelete)
  
  network4 <- collapseGraphs(graphs[7:8])
  toDelete <- which(igraph::degree(network4) == 0)
  network4 <- delete.vertices(network4, toDelete)
  # 
  #   network5 <- collapseGraphs(graphs[9:10])
  # toDelete <- which(igraph::degree(network5) == 0)
  #   network5 <- delete.vertices(network5, toDelete)
  # 
  #   
  #   print("----  time slice 1 similarities ----")
  minTime <- 1
  maxTime <- 14
  inactives <- which(igraph::degree(network1) == 0)
  userModels1 <- lapply(V(network1)$label, getUserModelHSHG, postData, network1, minTime, maxTime)
  #   #save(userModels1, file="user_model_dynamic1.Rdata")
  # 
  #userSemSimilarity1 <- getUserSimilarityMatrix(userModels1, threadSimMatrix, 0, 1, 0, inactives)
  #userSemSimilarity1[which(userSemSimilarity1 == Inf)] <- 0
  #save(userSemSimilarity1, file="user_sim_dynamic_finance_slice1.Rdata")
  load("user_sim_dynamic_finance_slice1.Rdata")
  userRegSimilarity1 <- REGE.for(get.adjacency(network1, sparse=FALSE))$E
  userStrSimilarity1 <- 1 - sedist(get.adjacency(network1, sparse=FALSE), method="euclidean")
  #   
  #   print("----  time slice 2 similarities ----")
  minTime <- 15
  maxTime <- 29
  inactives <- which(igraph::degree(network2) == 0)
  userModels2 <- lapply(V(network2)$label, getUserModelHSHG, postData, network2, minTime, maxTime)
  #   # save(userModels2, file="user_model_dynamic2.Rdata")
  #   
  #   userSemSimilarity2 <- getUserSimilarityMatrix(userModels2, threadSimMatrix, 0, 1, 0, inactives)
  #   userSemSimilarity2[which(userSemSimilarity2 == Inf)] <- 0  
  #   save(userSemSimilarity2, file="user_sim_dynamic_finance_slice2.Rdata")
  load("user_sim_dynamic_finance_slice2.Rdata")
  userRegSimilarity2 <- REGE.for(get.adjacency(network2, sparse=FALSE))$E
  userStrSimilarity2 <- 1 - sedist(get.adjacency(network2, sparse=FALSE), method="euclidean")
  #   
  #   print("----  time slice 3 similarities ----")
  minTime <- 30
  maxTime <- 44
  inactives <- which(igraph::degree(network3) == 0)
  userModels3 <- lapply(V(network3)$label, getUserModelHSHG, postData, network3, minTime, maxTime)
  #   #   save(userModels3, file="user_model_dynamic3.Rdata")
  #   
  #   userSemSimilarity3 <- getUserSimilarityMatrix(userModels3, threadSimMatrix, 0, 1, 0, inactives)
  #   userSemSimilarity3[which(userSemSimilarity3 == Inf)] <- 0
  #   save(userSemSimilarity3, file="user_sim_dynamic_finance_slice3.Rdata")
  load("user_sim_dynamic_finance_slice3.Rdata")
  userRegSimilarity3 <- REGE.for(get.adjacency(network3, sparse=FALSE))$E
  userStrSimilarity3 <- 1 - sedist(get.adjacency(network3, sparse=FALSE), method="euclidean")
  # 
  #   print("----  time slice 4 similarities ----")
  minTime <- 45
  maxTime <- 59
  inactives <- which(igraph::degree(network4) == 0)
  userModels4 <- lapply(V(network4)$label, getUserModelHSHG, postData, network4, minTime, maxTime)
  #   #   save(userModels4, file="user_model_dynamic4.Rdata")
  
  #   userSemSimilarity4 <- getUserSimilarityMatrix(userModels4, threadSimMatrix, 0, 1, 0)
  #   userSemSimilarity4[which(userSemSimilarity4 == Inf)] <- 0
  #   save(userSemSimilarity4, file="user_sim_dynamic_finance_slice4.Rdata")
  load("user_sim_dynamic_finance_slice4.Rdata")
  userRegSimilarity4 <- REGE.for(get.adjacency(network4, sparse=FALSE))$E
  userStrSimilarity4 <- 1 - sedist(get.adjacency(network4, sparse=FALSE), method="euclidean")
  
  
  #### Blockmodelling ###
  
  
  # Time slice 1
  
  simCor <- getSimilarityCorrelations(userStrSimilarity1, userRegSimilarity1, userSemSimilarity1)
  write.csv(simCor$r, paste("results_finance_001/slice1/similarity_correlations_", outputFilePrefix, ".csv", sep=""))
  write.csv(simCor$p, paste("results_finance_001/slice1/static_similarity_correlations_p_", outputFilePrefix, ".csv", sep=""))
  
  # --------------------------------------------------
  # semantic similarity = 0, regular similarity = 1
  # --------------------------------------------------
  
  bm <- getSemanticBlockmodel(network1, postData, threadSimMatrix, 0, 0, userSimMatrix=userRegSimilarity1)
  distances <- as.dist(1 - userSemSimilarity1)
  
  bmstats_static_s0_r1 <- evaluateBlockmodel(userModels1, keywords, bm, distances, paste("results_finance_001/slice1/s0_r1", outputFilePrefix, sep="_"))
  
  
  # --------------------------------------------------
  # semantic similarity = 1, regular similarity = 0
  # --------------------------------------------------
  
  bm <- getSemanticBlockmodel(network1, postData, threadSimMatrix, 0, 0, userSimMatrix=userSemSimilarity1)
  distances <- as.dist(1 - userSemSimilarity1)
  
  bmstats_static_s1_r0 <- evaluateBlockmodel(userModels1, keywords, bm, distances, paste("results_finance_001/slice1/s1_r0", outputFilePrefix, sep="_"))
  
  
  bmStatsTable1 <- rbind(bmstats_static_s0_r1, bmstats_static_s1_r0)
  
  row.names(bmStatsTable1) <- c("reg_only1", "sem_only1")
  
  
  # Time slice 2
  
  simCor <- getSimilarityCorrelations(userStrSimilarity2, userRegSimilarity2, userSemSimilarity2)
  write.csv(simCor$r, paste("results_finance_001/slice2/similarity_correlations_", outputFilePrefix, ".csv", sep=""))
  write.csv(simCor$p, paste("results_finance_001/slice2/static_similarity_correlations_p_", outputFilePrefix, ".csv", sep=""))
  
  # --------------------------------------------------
  # semantic similarity = 0, regular similarity = 1
  # --------------------------------------------------
  
  bm <- getSemanticBlockmodel(network2, postData, threadSimMatrix, 0, 0, userSimMatrix=userRegSimilarity2)
  distances <- as.dist(1 - userSemSimilarity2)
  
  bmstats_static_s0_r1 <- evaluateBlockmodel(userModels2, keywords, bm, distances, paste("results_finance_001/slice2/s0_r1", outputFilePrefix, sep="_"))
  
  
  # --------------------------------------------------
  # semantic similarity = 1, regular similarity = 0
  # --------------------------------------------------
  
  bm <- getSemanticBlockmodel(network2, postData, threadSimMatrix, 0, 0, userSimMatrix=userSemSimilarity2)
  distances <- as.dist(1 - userSemSimilarity2)
  
  bmstats_static_s1_r0 <- evaluateBlockmodel(userModels2, keywords, bm, distances, paste("results_finance_001/slice2/s1_r0", outputFilePrefix, sep="_"))
  bmStatsTable2 <- rbind(bmstats_static_s0_r1, bmstats_static_s1_r0)
  
  row.names(bmStatsTable2) <- c("reg_only2", "sem_only2")
  
  
  # Time slice 3
  
  simCor <- getSimilarityCorrelations(userStrSimilarity3, userRegSimilarity3, userSemSimilarity3)
  write.csv(simCor$r, paste("results_finance_001/slice3/similarity_correlations_", outputFilePrefix, ".csv", sep=""))
  write.csv(simCor$p, paste("results_finance_001/slice3/static_similarity_correlations_p_", outputFilePrefix, ".csv", sep=""))
  
  # --------------------------------------------------
  # semantic similarity = 0, regular similarity = 1
  # --------------------------------------------------
  
  bm <- getSemanticBlockmodel(network3, postData, threadSimMatrix, 0, 0, userSimMatrix=userRegSimilarity3)
  distances <- as.dist(1 - userSemSimilarity3)
  
  bmstats_static_s0_r1 <- evaluateBlockmodel(userModels3, keywords, bm, distances, paste("results_finance_001/slice3/s0_r1", outputFilePrefix, sep="_"))
  
  
  # --------------------------------------------------
  # semantic similarity = 1, regular similarity = 0
  # --------------------------------------------------
  
  bm <- getSemanticBlockmodel(network3, postData, threadSimMatrix, 0, 0, userSimMatrix=userSemSimilarity3)
  distances <- as.dist(1 - userSemSimilarity3)
  
  bmstats_static_s1_r0 <- evaluateBlockmodel(userModels3, keywords, bm, distances, paste("results_finance_001/slice3/s1_r0", outputFilePrefix, sep="_"))
  bmStatsTable3 <- rbind(bmstats_static_s0_r1, bmstats_static_s1_r0)
  
  row.names(bmStatsTable3) <- c("reg_only3", "sem_only3")
  
  
  # Time slice 4
  
  simCor <- getSimilarityCorrelations(userStrSimilarity4, userRegSimilarity4, userSemSimilarity4)
  write.csv(simCor$r, paste("results_finance_001/slice4/similarity_correlations_", outputFilePrefix, ".csv", sep=""))
  write.csv(simCor$p, paste("results_finance_001/slice4/static_similarity_correlations_p_", outputFilePrefix, ".csv", sep=""))
  
  # --------------------------------------------------
  # semantic similarity = 0, regular similarity = 1
  # --------------------------------------------------
  
  bm <- getSemanticBlockmodel(network4, postData, threadSimMatrix, 0, 0, userSimMatrix=userRegSimilarity4)
  distances <- as.dist(1 - userSemSimilarity4)
  
  bmstats_static_s0_r1 <- evaluateBlockmodel(userModels4, keywords, bm, distances, paste("results_finance_001/slice4/s0_r1", outputFilePrefix, sep="_"))
  
  
  # --------------------------------------------------
  # semantic similarity = 1, regular similarity = 0
  # --------------------------------------------------
  
  bm <- getSemanticBlockmodel(network4, postData, threadSimMatrix, 0, 0, userSimMatrix=userSemSimilarity4)
  distances <- as.dist(1 - userSemSimilarity4)
  
  bmstats_static_s1_r0 <- evaluateBlockmodel(userModels4, keywords, bm, distances, paste("results_finance_001/slice4/s1_r0", outputFilePrefix, sep="_"))
  
  bmStatsTable4 <- rbind(bmstats_static_s0_r1, bmstats_static_s1_r0)
  
  row.names(bmStatsTable4) <- c("reg_only4", "sem_only4")
  
  
  dynamicStats <- rbind(bmStatsTable1, bmStatsTable2, bmStatsTable3, bmStatsTable4)
  write.csv(dynamicStats, "results_finance_001/dynamics.csv")
  
  
  
}

evaluateDynamics <- function() {
  
  postData <- read.csv2("../../network_extraction/post_tables/finance_001_featurised/latest/classified.csv")
  keywords <- read.csv("../finance_001/thread_keywords_finance_001/calais_annotations/calais_tags_finance.csv", as.is="keyword")
  threadSimMatrix <- readMM("../finance_001/thread_sem_matrix_wordnet.txt")
  dimNames <- read.table("../finance_001/thread_sem_sim_wordnet_dimnames.txt", as.is = "x")$x
  rownames(threadSimMatrix) <- dimNames
  colnames(threadSimMatrix) <- dimNames
  threadSimMatrix[which(threadSimMatrix > 1)] <- threadSimMatrix[which(threadSimMatrix > 1)] / 1e+308
  
  #   threadSimMatrix <- readMM("../finance_001/thread_sem_matrix_lsa.txt")
  #   dimNames <- read.table("../finance_001/thread_sem_sim_lsa_dimnames.txt", as.is = "x")$x
  #   rownames(threadSimMatrix) <- dimNames
  #   colnames(threadSimMatrix) <- dimNames
  
  graphDir <- "../../network_extraction/classified_datasets/fincance_001/latest/timeslices/"
  graphFiles <- list.files(graphDir, full.names = TRUE)
  graphs <- lapply(graphFiles, function(file) {
    
    read.graph(file, "gml")
  })
  
  network <- collapseGraphs(graphs)
  toDelete <- which(igraph::degree(network) == 0)
  network <- delete.vertices(network, toDelete)
  
  outputFilePrefix <- "finance_001"
  
  load("results_finance_001/user_sim_static.Rdata")
  userSemSim[which(userSemSim == Inf)] <- 0
  userRegSim <- REGE.for(get.adjacency(network, sparse=FALSE))$E
  
  network1 <- collapseGraphs(graphs[1:2])
  
  network2 <- collapseGraphs(graphs[3:4])
  
  network3 <- collapseGraphs(graphs[5:6])
  
  network4 <- collapseGraphs(graphs[7:8])
  
  
  minTime <- 1
  maxTime <- 14
  userModels1 <- lapply(V(network1)$label, getUserModelHSHG, postData, network1, minTime, maxTime)
  
  minTime <- 15
  maxTime <- 29
  userModels2 <- lapply(V(network2)$label, getUserModelHSHG, postData, network2, minTime, maxTime)
  
  minTime <- 30
  maxTime <- 44
  userModels3 <- lapply(V(network3)$label, getUserModelHSHG, postData, network3, minTime, maxTime)
  
  minTime <- 45
  maxTime <- 59
  userModels4 <- lapply(V(network4)$label, getUserModelHSHG, postData, network4, minTime, maxTime)
  
  # ---------------  
  # dynamics
  # ---------------
  userSimMatrix <- (2 * userSemSim + userRegSim) / 3
  distances <- as.dist(1 - userSimMatrix)
  h <- hclust(distances, method = "ward.D")
  
  numClusters <- nselectboot(distances, clustermethod=disthclustCBI, method="ward.D")$kopt
  
  part <- cutree(h, numClusters)
  
  # isolated nodes should not contribute to the denstity
  nodesConnected1 <- length(which(igraph::degree(network1) > 0))
  dens1 <- length(E(network1)) / (nodesConnected1 * (nodesConnected1 - 1))
  nodesConnected2 <- length(which(igraph::degree(network2) > 0))
  dens2 <- length(E(network2)) / (nodesConnected2 * (nodesConnected2 - 1))
  nodesConnected3 <- length(which(igraph::degree(network3) > 0))
  dens3 <- length(E(network3)) / (nodesConnected3 * (nodesConnected3 - 1))
  nodesConnected4 <- length(which(igraph::degree(network4) > 0))
  dens4 <- length(E(network4)) / (nodesConnected4 * (nodesConnected4 - 1))
  
  bm1 <- crit.fun(get.adjacency(network1, sparse=FALSE), part, approach="bin", blocks=c("null","reg"), blockWeights=c(null=1,reg=dens1/(1-dens1)), norm=TRUE)
  bm2 <- crit.fun(get.adjacency(network2, sparse=FALSE), part, approach="bin", blocks=c("null","reg"), blockWeights=c(null=1,reg=dens2/(1-dens2)), norm=TRUE)
  bm3 <- crit.fun(get.adjacency(network3, sparse=FALSE), part, approach="bin", blocks=c("null","reg"), blockWeights=c(null=1,reg=dens3/(1-dens3)), norm=TRUE)
  bm4 <- crit.fun(get.adjacency(network4, sparse=FALSE), part, approach="bin", blocks=c("null","reg"), blockWeights=c(null=1,reg=dens4/(1-dens4)), norm=TRUE)
  
  bmstats_dynamic_s1_r2 <- evaluateDynamicBlockmodel(list(userModels1, userModels2, userModels3, userModels4), keywords, list(bm1, bm2, bm3, bm4), distances, 
                                                     list(network1, network2, network3, network4), paste("results_finance_001/dynamics/dynamic_s2_r1", outputFilePrefix, sep="_"))  
}

evaluateRandomPartition <- function() {
  
  postData <- read.csv2("../../network_extraction/post_tables/finance_001_featurised/latest/classified.csv")
  keywords <- read.csv("../finance_001/thread_keywords_finance_001/calais_annotations/calais_tags_finance.csv", as.is="keyword")
  threadSimMatrix <- readMM("../finance_001/thread_sem_matrix_wordnet.txt")
  dimNames <- read.table("../finance_001/thread_sem_sim_wordnet_dimnames.txt", as.is = "x")$x
  rownames(threadSimMatrix) <- dimNames
  colnames(threadSimMatrix) <- dimNames
  threadSimMatrix[which(threadSimMatrix > 1)] <- threadSimMatrix[which(threadSimMatrix > 1)] / 1e+308
  
  #   threadSimMatrix <- readMM("../finance_001/thread_sem_matrix_lsa.txt")
  #   dimNames <- read.table("../finance_001/thread_sem_sim_lsa_dimnames.txt", as.is = "x")$x
  #   rownames(threadSimMatrix) <- dimNames
  #   colnames(threadSimMatrix) <- dimNames
  
  graphDir <- "../../network_extraction/classified_datasets/fincance_001/latest/timeslices/"
  graphFiles <- list.files(graphDir, full.names = TRUE)
  graphs <- lapply(graphFiles, function(file) {
    
    read.graph(file, "gml")
  })
  
  network <- collapseGraphs(graphs)
  toDelete <- which(igraph::degree(network) == 0)
  network <- delete.vertices(network, toDelete)
  
  outputFilePrefix <- "finance_001"
  
  ####################################################
  #
  # Static blockmodels
  #
  ####################################################
  
  minTime <- 1
  maxTime <- 70
  print("user models static")
  userModels <- lapply(V(network)$label, getUserModelHSHG, postData, network, minTime, maxTime)
  #     #save(userModels, file="user_model_static.Rdata")
  #     print("user similarity static")
  #     #userSemSim <- getUserSimilarityMatrix(userModels, threadSimMatrix, 0, 1, 0)
  #     save(userSemSim, file="results_finance_001/user_sim_static.Rdata")
  load("results_finance_001/user_sim_static.Rdata")
  userSemSim[which(userSemSim == Inf)] <- 0
  userRegSim <- REGE.for(get.adjacency(network, sparse=FALSE))$E
  userStrSim <- 1 - sedist(get.adjacency(network, sparse=FALSE), method="euclidean")
  
  nodesConnected <- length(which(igraph::degree(network) > 0))
  dens <- length(E(network)) / (nodesConnected * (nodesConnected - 1))
  distances <- as.dist(1 - userSemSim)
  sapply(1:10, function(i) {
    
    part <- genRandomPar(3, length(V(network)), seed = i)
    bm <- crit.fun(get.adjacency(network, sparse=FALSE), part, approach="bin", blocks=c("null","reg"), blockWeights=c(null=1,reg=dens/(1-dens)), norm=TRUE)
    evaluateBlockmodel(userModels, keywords, bm, distances, paste("results_finance_001/static_s0_r1", outputFilePrefix, sep="_"))
  })
}