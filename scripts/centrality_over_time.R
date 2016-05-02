require(igraph)
require(dtw)
require(fpc)
require(cluster)
require(psych)
require(ggplot2)

occurencesOverTime <- function(timeslices) {
  
  rowSums(sapply(timeslices, function(ts) {
    
    (degree(ts) > 0)
  }))
}

usersOverTime <- function(timeslices) {
  
  colSums(sapply(timeslices, function(ts) {
    
    (degree(ts) > 0)
  }))
}

centralityOverTime <- function(timeslices, centralityFunction=degree, dumpingFactor=1, dumpingFunction="linear") {
  
  sapply(1:length(timeslices), function(i) {
    
    graph <- timeslices[[i]]
    if (is.null(E(graph)$weight)) {
      
      wAdj <- get.adjacency(graph)
    } else {
      
      wAdj <- get.adjacency(graph, attr = "weight")  
    }
    
    if (i > 1) {
      
      sapply(1:(i-1), function(j) {
      
        if (!is.null(dumpingFactor)) {
          
          if (dumpingFunction == "linear") {
            
            wAdj <<- wAdj + (get.adjacency(timeslices[[j]]) * dumpingFactor^(i-j))  
      
          } else if (dumpingFunction == "inv_exp") {
            
            wAdj <<- wAdj + exp(j - i) * get.adjacency(timeslices[[j]])  
          }
          
        } else {
          
          wAdj <<- wAdj + get.adjacency(timeslices[[j]])
        } 
      })
    }
    gAggregated <- graph.adjacency(wAdj, weighted = TRUE)
    
    centralityFunction(gAggregated)
  })
}

outDegree <- function(graph) {
  
  graph.strength(graph, mode="out")
}

inDegree <- function(graph) {
  
  graph.strength(graph, mode="in")
}

expertiseRank <- function(graph) {
  
  
}

outreach <- function(graph, vertices=V(graph)) {
  
  numAnswers <- graph.strength(graph, mode="out")
  
  sapply(vertices, function(v) {
    
    incEdges <- incident(graph, v, mode="out")
    
    ent <- 0
    if ((length(incEdges) > 0) && (numAnswers[v] > 1)) {
      sapply(incEdges$weight, function(w) {
        
        p <- w / numAnswers[v]
        ent <<- ent - (p * log(p))
      })
    }
    numAnswers[v] * (ent + 1)
  })
}

inreach <- function(graph, vertices=V(graph)) {
  
  numAnswers <- graph.strength(graph, mode="in")
  
  sapply(vertices, function(v) {
    
    incEdges <- incident(graph, v, mode="in")
    
    ent <- 0
    if ((length(incEdges) > 0) && (numAnswers[v] > 1)) {
      sapply(incEdges$weight, function(w) {
        
        p <- w / numAnswers[v]
        ent <<- ent - (p * log(p))
      })
    }
    
    numAnswers[v] * (ent + 1)
  })
}

zScore <- function(graph) {
  
  numQuestions <- graph.strength(graph, mode="in")
  numAnswers <- graph.strength(graph, mode="out")
  
  sapply(V(graph), function(v) {
    
    (numAnswers[v] - numQuestions[v]) / sqrt((numAnswers[v] + numQuestions[v]))
  })
}

hubScore <- function(graph) {
  
  hub.score(graph)$vector
}

authorityScore <- function(graph) {
  
  authority.score(graph)$vector
}

collapseGraphs <- function(graphs) {
  
  isFirst <- TRUE
  adj <- NULL
  
  sapply(1:length(graphs), function(i) {
    
    if (isFirst) {
      
      isFirst <<- FALSE
      adj <<- get.adjacency(graphs[[i]])
    } else {
      
      adj <<- adj + get.adjacency(graphs[[i]])
    }
  })
  
  g <- graph.adjacency(adj, weighted = TRUE)
  V(g)$label <- V(graphs[[1]])$label
  g
}

startCentralityOverTime <- function() {
  
  graphDir <- "data/timeslices"
  
  graphFiles <- list.files(graphDir, full.names = TRUE)
  
  graphs <- lapply(graphFiles, read.graph, "gml")
  
  names <- V(graphs[[1]])$label
  
  #graphs <- list(collapseGraphs(graphs)) #####
  
  odCentrality <- centralityOverTime(graphs, outDegree, dumpingFactor=0.9, dumpingFunction = "linear")
  row.names(odCentrality) <- names
  write.csv(odCentrality, "outdegree_over_time.csv")
  idCentrality <- centralityOverTime(graphs, inDegree, dumpingFactor=0.9, dumpingFunction = "linear")
  row.names(idCentrality) <- names
  write.csv(idCentrality, "indegree_over_time.csv")
  btwCentrality <- centralityOverTime(graphs, betweenness, dumpingFactor=0.9, dumpingFunction = "linear")
  row.names(btwCentrality) <- names
  write.csv(btwCentrality, "betweenness_over_time.csv")
  outreach <- centralityOverTime(graphs, outreach, dumpingFactor=0.9, dumpingFunction = "linear")
  row.names(outreach) <- names
  write.csv(outreach, "outreach_over_time.csv")
  inreach <- centralityOverTime(graphs, inreach, dumpingFactor=0.9, dumpingFunction = "linear")
  row.names(inreach) <- names
  write.csv(inreach, "inreach_over_time.csv")
  z <- centralityOverTime(graphs, zScore, dumpingFactor=0.9, dumpingFunction = "linear")
  row.names(z) <- names
  write.csv(z, "z-score_over_time.csv")
  hubs <- centralityOverTime(graphs, hubScore, dumpingFactor=0.9, dumpingFunction = "linear")
  row.names(hubs) <- names
  write.csv(hubs, "hub-score_over_time.csv")
  authorities <- centralityOverTime(graphs, authorityScore, dumpingFactor=0.9, dumpingFunction = "linear")
  row.names(authorities) <- names
  write.csv(authorities, "authorities_over_time.csv")
  
  list(od=odCentrality, id=idCentrality, btw=btwCentrality, inreach=inreach, outreach=outreach, hub_score=hubs, authority_score=authrities, z_score=z)
}

plotInAndOutreachTrajectories <- function(inreachData, outreachData, medoids) {
  
  inreachCol <- c()
  outreachCol <- c()
  clustNames <- c()
  times <- c()
  i <- 1
  sapply(medoids, function(medoid) {
    
    row <- which(rownames(inreachData) == medoid)
    inreachCol <<- append(inreachCol, as.numeric(inreachData[row,]))
    outreachCol <<- append(outreachCol, as.numeric(outreachData[row,]))
    clustNames <<- append(clustNames, paste("cf", rep(i, ncol(inreachData))))
    times <<- append(times, 1:ncol(inreachData))
    i <<- i + 1
  })
  values <- append(inreachCol, outreachCol)
  measures <- append(rep("inreach", length(inreachCol)),rep("outreach", length(inreachCol)))
  clusters <- append(clustNames, clustNames)
  slices <- append(times, times)                   
  dataTable <- data.frame(value=values, time_slice=slices, measure=measures, group=clusters)
  
  diagram <- ggplot(data=dataTable, aes(x=time_slice, y=value, group = measure, colour = measure, shape=measure, fill="white")) +
    geom_line() +
    geom_point( size=4, fill="white") + facet_wrap(~group, ncol = 2, scales = "free") + theme(axis.text=element_text(size=18)) + theme_bw()
  
  plot(diagram)
  
  diagram
}

# This function clusters users according to their in- and outreach trajectories over time and plots the result.
getInAndOutreachTrajectories <- function(inreachData, outreachData, numClusters=15, iterations=15) {
  
  # clean data
  toKeep <- union(which(rowSums(inreachData) > 0), which(rowSums(outreachData) > 0))
  print(paste("Keep", length(toKeep), "of", nrow(inreachData), "examples"))
  inreachData <- inreachData[toKeep,]
  outreachData <- outreachData[toKeep,]
  
  
  inrDist <- dist(inreachData, method = "manhattan")
  outrDist <- dist(outreachData, method = "manhattan")
  cbdt <- as.dist((inrDist + outrDist) / 2)
  

  # Alternating k-medoids
  print("-- Alternating k-medoids")
  i <- 1
  clustering1 <- pam(inrDist, k = numClusters, diss = TRUE)
  
  
  while(i < iterations) {
    
    initMedoids <- which(rownames(outreachData) %in% clustering1$medoids)
    clustering2 <- pam(outrDist, k = numClusters, medoids = initMedoids, diss = TRUE)
    initMedoids <- which(rownames(inreachData) %in% clustering2$medoids)
    clustering1 <- pam(inrDist, k = numClusters, medoids = initMedoids, diss = TRUE)
    
    cstats <- cluster.stats(clustering = clustering1$clustering, alt.clustering = clustering2$clustering, compareonly = TRUE)
    print(paste("Run", i, "Clustering similarity Rand index =", cstats$corrected.rand, "vi =", cstats$vi))
    i <- i + 1
  }
  altKmedoidsResult <- clustering1
  inStats <- cluster.stats(inrDist, clustering = altKmedoidsResult$clustering)
  outStats <- cluster.stats(outrDist, clustering = altKmedoidsResult$clustering)
  print(paste("in dunn:", inStats$dunn2, "out dunn:", outStats$dunn2, "sum:", (inStats$dunn2 + outStats$dunn2)))
  print(paste("in bw:", (inStats$average.between / inStats$average.within), "out bw:", outStats$average.between / outStats$average.within, 
              "sum:", ((inStats$average.between / inStats$average.within) + (outStats$average.between / outStats$average.within))))
  
  plotInAndOutreachTrajectories(inreachData, outreachData, clustering1$medoids) 

  
  altKmedoidsResult
}

# Starts plotting inreach and outreach of users over time.
startTrajectoryPlot <- function() {
    
  inreachData <- read.csv("data/inreach_over_time.csv", row.names=1)
  outreachData <- read.csv("data/outreach_over_time.csv", row.names=1)
  
  getInAndOutreachTrajectories(inreachData, outreachData)
}
