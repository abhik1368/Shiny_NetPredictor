library(bipartite)
data(Safariland)
plotweb(Safariland)
visweb(Safariland)
networklevel(Safariland)
specieslevel(Safariland)
closeness_w(Safariland)
web <- matrix(0, 10,10)
web[1,1:3] <- 1
web[2,4:5] <- 1
web[3:7, 6:8] <- 1
web[8:10, 9:10] <- 1
web <- web[-c(4:5),] #oh, and make it asymmetric!
web <- web[,c(1:5, 9,10, 6:8)] #oh, and make it non-diagonal

library(visNetwork)
nodes <- jsonlite::fromJSON("https://raw.githubusercontent.com/DataKnowledge/DataKnowledge.github.io/master/visNetwork/data/nodes_miserables.json")

edges <- jsonlite::fromJSON("https://raw.githubusercontent.com/DataKnowledge/DataKnowledge.github.io/master/visNetwork/data/edges_miserables.json")


visNetwork(nodes, edges, height = "700px", width = "100%") %>%
    visOptions(selectedBy = "group", highlightNearest = TRUE, nodesIdSelection = TRUE)


mynet <- read.csv("Results-3.txt",header=TRUE,sep='\t')
nodes <- unique(c(as.character(mynet$Proteins), as.character(mynet$Drugs)))
group1 <- length(unique(mynet$Proteins))
group2 <- length(unique(mynet$Drugs))
id=seq(0,length(nodes)-1,1)
label=nodes

nodeData <- data.frame(id ,label,group=c(rep("Proteins",group1),rep("Drugs",group2)),stringsAsFactors=FALSE)
 

#nodeData <- data.frame(name,stringsAsFactors=FALSE)
#nodeData <- data.frame(name,group=c(rep(10,group1),rep(25,group2)),stringsAsFactors=FALSE)
colnames(mynet)[1] <- "from"
colnames(mynet)[2] <- "to"
edgeList <- mynet[, c("from","to")]
edgeList$dashes <- ifelse(mynet$type == "True Interactions",FALSE,TRUE)
#nodesID <- data.frame(nodes,id=seq(0,length(nodes)-1,1))

edgeList$from <- with(nodeData, id[match(edgeList$from, nodes)])
edgeList$to <- with(nodeData, id[match(edgeList$to,nodes)])

visNetwork(nodeData, edgeList,height="700px",width = "100%") %>% visInteraction(dragNodes = FALSE,hover = TRUE,navigationButtons = TRUE,tooltipDelay = 0) %>% visLegend() %>% visNodes(scaling=list(min=20),font=list(size=24)) %>%
visOptions(selectedBy = "group", highlightNearest = TRUE, nodesIdSelection = TRUE) %>% visPhysics(stabilization=FALSE) %>% visLayout(randomSeed = 123) %>% visPhysics(stabilization=list(iterations=200))


edgeList$dashes <- ifelse(mynet$type == "True Interactions","FALSE","TRUE")
    