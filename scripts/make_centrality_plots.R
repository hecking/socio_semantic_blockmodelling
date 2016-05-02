require(igraph)
source("centrality_over_time.R")

graphDir <- "data/timeslices/"
graphFiles <- list.files(graphDir)

graphs <- lapply(graphFiles, function(file) {
  
  read.graph(paste(graphDir, file, sep=""), "gml")
})

g <- read.graph("fincance_001/latest/qa_net_finance_001.gml", "gml")  

odCentrality <- centralityOverTime(graphs, outDegree, dumpingFactor=0.95, dumpingFunction = "linear")
idCentrality <- centralityOverTime(graphs, inDegree, dumpingFactor=0.95, dumpingFunction = "linear")
btwCentrality <- centralityOverTime(graphs, betweenness, dumpingFactor=0.95, dumpingFunction = "linear")
outreach <- centralityOverTime(graphs, outreach, dumpingFactor=0.95, dumpingFunction = "linear")
inreach <- centralityOverTime(graphs, inreach, dumpingFactor=0.95, dumpingFunction = "linear")

hubs <- which(degree(g) > 15)

out <- NULL
hub <- NULL
sapply(hubs, function(node) {
  par(mar=c(3, 4, 2, 2), xpd=TRUE)
  plot(inreach[node,], type="o", col="blue", ann=FALSE, ylim=c(0, (max(max(inreach[node,]), max(outreach[node,])) + 5)))
  
  box()
  
  lines(outreach[node,], type="o", pch=22, lty=2, col="red")
  
  
  #title(main="In- and Out-reach Over Time", font.main=4)
  
  title(xlab="Time slice", col.lab=rgb(0,0,0), line=1.6, cex.lab=1.2)
  title(ylab="In-/Outreach", col.lab=rgb(0,0,0), cex.lab=1.2, line=1.8)
}) 