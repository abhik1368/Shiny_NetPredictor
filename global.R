# Copyright (C) 2015 Abhik Seal <abhik1368@gmail.com>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#################################################################################
## INSTALLATION OF PACKAGES FIRST
# list.of.packages <- c("devtools","DBI","RSQLite","shinyBS","doParallel","foreach","netpredictor","lpbrim","shinysky","shinyjs","shinythemes","reshape2","rlist","htmltools","igraph","gdata","shiny","data.table","visNetwork","DBI","RSQLite") # replace xx and yy with package names
#
#
#  packages.auto <- function(x) {
#     x <- as.character(substitute(x))
#     if(isTRUE(x %in% .packages(all.available=TRUE))) {
#         eval(parse(text = sprintf("require(\"%s\")", x)))
#      } else {
#          eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
#     }
#      if(isTRUE(x %in% .packages(all.available=TRUE))) {
#         eval(parse(text = sprintf("require(\"%s\")", x)))
#      } else {
#
#         devtools::install_github("abhik1368/netpredictor")
#         devtools::install_github("AnalytixWare/ShinySky")
#          devtools::install_github("PoisotLab/lpbrim")
#          devtools::install_github("shiny-gridster", "wch")
#          devtools::install_github("ebailey78/shinyBS")
#
#
#     }
# }
#  bio.packages <- c("biomaRt","ReactomePA","clusterProfiler")
#
#  bioconduct.auto <- function(x) {
#     source("https://bioconductor.org/biocLite.R")
#     x <- as.character(substitute(x))
#     if(isTRUE(x %in% .packages(all.available=TRUE))) {
#          eval(parse(text = sprintf("require(\"%s\")", x)))
#      } else {
#         #update.packages(ask= FALSE) #update installed packages.
#         eval(parse(text = sprintf("biocLite(\"%s\", dependencies = TRUE)", x)))
#      }
#  }
#  packages.auto(list.of.packages)
#  bioconduct.auto(bio.packages)
#################################################################################

#################################################################################
## Load the libraries
library(shiny)
library(dplyr)
library(RSQLite)
library(netpredictor)
library(igraph)
library(shiny)
require(shinysky)
library(reshape2)
library(lpbrim)
library(rlist)
library(doParallel)
library(foreach)
library(org.Hs.eg.db)
library(tidyr)
library(DBI)
library(biomaRt)
set.seed(12345)
#################################################################################
# Libraries ---------------------------------------------------------------
# library("sodium")
# library("dtplyr")
library(shinyWidgets)
# Global Variables --------------------------------------------------------

# app_name <- "User Auth"
# database <- "authorization.db"
# support_contact <- "abhik1368@gmail.com"

## db connection
con2 <- dbConnect(SQLite(), "data/ppi.sqlite")
Q1 <- sprintf("SELECT * FROM human_ppi")
interactome_data <- dbGetQuery(con2,Q1)
dbDisconnect(con2)
human <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

load("data/UnWeighted_ConsensusPath.rda")
load("data/Weighted_ConsensusPath.rda")
load("data/Weighted_String.rda")
load("data/UnWeighted_String.rda")
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#################################################################################
## Function for Netpredictor
#################################################################################

getProp <- function(data){

    ## degree Centrality of the Bipartite Graph
    deg <- get.biDegreeCentrality(data,SM=TRUE)

    ## Compute Graph density of Bipartite Graph
    den <- get.biDensity(data)

    ## Compute betweeness centrality of Bipartite Graph
    bwt <- get.biBetweenessCentrality(data)

    edges_count <- sum(data)

    no_drugs <- dim(data)[1]
    no_proteins <- dim(data)[2]

    NetworkProperties=c("Degree Centrality Protein","Degree Centrality Drugs",
                        "Density of Network" ,
                        "Betweeeness Centrality Protein",
                        "Betweeness Centrality Drugs",
                        "Number of Interactions",
                        "No of Drugs",
                        "No of proteins")
    Values=c(as.numeric(as.character(deg$SM.bi.degreeCentrality.columns[[1]])),
             as.numeric(as.character(deg$SM.bi.degreeCentrality.rows[[1]])),
             as.numeric(as.character(den[[1]])),
             as.numeric(as.character(bwt$SM.bi.BetweenessCentrality.columns[[1]])),
             bwt$SM.bi.BetweenessCentrality.rows[[1]],
             edges_count[[1]],
             as.numeric(as.character(dim(data)[1])),
             as.numeric(as.character(dim(data)[2])))
    df <- data.frame(NetworkProperties,Values)

        return (df)
#     bg=graph.incidence(A)
#     pr=bipartite.projection(bg)
#     gp <- get.Communities(pr$proj1)
}


getMod <- function(dt){
    M <- as.matrix(dt)
    M <- M[rowSums(M)>0, colSums(M)>0]
    mod <- findModules(M, iter=10, sparse=FALSE)
    d <- getmodules(mod)
    p <- list.filter(d, x ~ dim(x)[1] > 0  )
    p

}


getheat <- function(data,Alpha,Lambda){

    load(data)

    A <- t(adjm)
    R <- nbiNet(A,alpha=Alpha,lamda=Lambda,format="matrix")
    tab <- getTopresults(A,R,top=10,druglist=NULL)
    colnames(tab)[1] <- "Proteins"
    colnames(tab)[2] <- "Drugs"
    tab$Proteins <- as.character(tab$Proteins)
    tab$Drugs <- as.character(tab$Drugs)
    tab$score <- as.numeric(as.character(tab$score))
    res <- tab[tab$score > 0, ]

    return(res)

    }

getnbi <- function(data,Alpha,Lambda){

      load(data)
      A <- t(adjm)
      S1 = as.matrix(compSim)
      S2 = as.matrix(protSim)


      # Use min/max and setProgress(value=x) for progress bar

      P1 <- nbiNet(A,alpha=Alpha, lamda=Lambda,  s1=S1, s2=S2,format = "matrix")
      tab <- getTopresults(A,P1,top=10,druglist=NULL)
      colnames(tab)[1] <- "Proteins"
      colnames(tab)[2] <- "Drugs"
      tab$Proteins <- as.character(tab$Proteins)
      tab$Drugs <- as.character(tab$Drugs)
      tab$score <- as.numeric(as.character(tab$score))
      res <- tab[tab$score > 0, ]

      return(res)

}

getrwr <- function(data,c){

    load(data)
    A <- t(adjm)
    S1 = as.matrix(compSim)
    S2 = as.matrix(protSim)
    g1 = graph.incidence(A)
    Q = biNetwalk(g1,s1=S1,s2=S2,normalise="laplace",dataSeed=NULL,restart=c,verbose=T)
    tab <- getTopresults(A,Q,top=10,druglist=NULL)
    colnames(tab)[1] <- "Proteins"
    colnames(tab)[2] <- "Drugs"
    tab$Proteins <- as.character(tab$Proteins)
    tab$Drugs <- as.character(tab$Drugs)
    tab$score <- as.numeric(as.character(tab$score))
    res <- tab[tab$score > 0, ]

    return(res)
}

getcombo <- function(data,Alpha,Lambda,C){

    load(data)
    A <- t(adjm)
    S1 = as.matrix(compSim)
    S2 = as.matrix(protSim)
    g1 = graph.incidence(A)
    busyIndicator("Running RWR",wait = 0)
    mat1 = biNetwalk(g1,s1=S1,s2=S2,normalise="laplace",dataSeed=NULL,restart=C,verbose=T)
    busyIndicator("Running NBI",wait = 0)
    mat2 <- nbiNet(A, lamda=Lambda, alpha=Alpha, s1=S1, s2=S2,format = "matrix")
    mat = (mat1+mat2)/2
    tab <- getTopresults(A,mat,top=10,druglist=NULL)
    colnames(tab)[1] <- "Proteins"
    colnames(tab)[2] <- "Drugs"
    tab$Proteins <- as.character(tab$Proteins)
    tab$Drugs <- as.character(tab$Drugs)
    tab$score <- as.numeric(as.character(tab$score))
    res <- tab[tab$score > 0, ]
    return(res)
}


getCustomHeat <- function(dataDT,heat_alpha,heat_lamda){

    A <- t(dataDT)

    R <- nbiNet(A,alpha=heat_alpha,lamda=heat_lamda,format="matrix")
    tab <- getTopresults(A,R,top=10,druglist=NULL)
    colnames(tab)[1] <- "Proteins"
    colnames(tab)[2] <- "Drugs"
    tab$Proteins <- as.character(tab$Proteins)
    tab$Drugs <- as.character(tab$Drugs)
    tab$score <- as.numeric(as.character(tab$score))
    res <- tab[tab$score > 0, ]

    return(res)
}

getCustomNBI <- function (dataDT,dataD,dataT,nbi_alpha,nbi_lamda){

    A <- t(dataDT)
    S1 = as.matrix(dataD)
    S2 = as.matrix(dataT)


    # Use min/max and setProgress(value=x) for progress bar

    P1 <- nbiNet(A,alpha=nbi_alpha, lamda=nbi_lamda,  s1=S1, s2=S2,format = "matrix")
    tab <- getTopresults(A,P1,top=10,druglist=NULL)
    colnames(tab)[1] <- "Proteins"
    colnames(tab)[2] <- "Drugs"
    tab$Proteins <- as.character(tab$Proteins)
    tab$Drugs <- as.character(tab$Drugs)
    tab$score <- as.numeric(as.character(tab$score))
    res <- tab[tab$score > 0, ]

    return(res)
}


getCustomRWR <- function (dataDT,dataD,dataT,C){

    A <- t(dataDT)
    S1 = as.matrix(dataD)
    S2 = as.matrix(dataT)
    g1 = graph.incidence(A)
    Q = biNetwalk(g1,s1=S1,s2=S2,normalise="laplace",dataSeed=NULL,restart=C, verbose=T)
    tab <- getTopresults(A,Q,top=10,druglist=NULL)
    colnames(tab)[1] <- "Proteins"
    colnames(tab)[2] <- "Drugs"
    tab$Proteins <- as.character(tab$Proteins)
    tab$Drugs <- as.character(tab$Drugs)
    tab$score <- as.numeric(as.character(tab$score))
    res <- tab[tab$score > 0, ]

    return(res)
}


getCustomCombo <- function (dataDT,dataD,dataT,nc_alpha,nc_lamda,C){

    A <- t(dataDT)
    S1 = as.matrix(dataD)
    S2 = as.matrix(dataT)
    g1 = graph.incidence(A)
    mat1 = biNetwalk(g1,s1=S1,s2=S2,normalise="laplace",dataSeed=NULL,restart=C,verbose=T)
    mat2 <- nbiNet(A, lamda=nc_lamda, alpha=nc_alpha, s1=S1, s2=S2,format = "matrix")
    mat = (mat1+mat2)/2
    tab <- getTopresults(A,mat,top=10,druglist=NULL)
    colnames(tab)[1] <- "Proteins"
    colnames(tab)[2] <- "Drugs"
    tab$Proteins <- as.character(tab$Proteins)
    tab$Drugs <- as.character(tab$Drugs)
    tab$score <- as.numeric(as.character(tab$score))
    res <- tab[tab$score > 0, ]
    return(res)
}

extractUC <- function(x.1,x.2,...){
    x.1p <- do.call("paste", x.1)
    x.2p <- do.call("paste", x.2)
    x.1[! x.1p %in% x.2p, ]
}

getPermScores <- function(sigNetwork,A,sig){
    sigNetwork[sigNetwork >= sig] <- 10
    PD <- melt(sigNetwork)
    P <- PD[ order(PD[,3]), ]
    FPD <- P[P$value!=10,]
    FP <- FPD[,1:2]
    ## Original dataframe
    OD <- melt(A)
    OD <- OD[OD$value == 1,]
    OD <- OD[,1:2]
    m <- extractUC(FP,OD)
    new <- merge(m,FPD)
    old <- merge(OD,FPD)
    new$outcome <- "Predicted Interactions"
    old$outcome <- "True Interactions"
    permResult <- rbind(new,old)
    return(permResult)
}

permTest <- function(dataDT,dataD,dataT,restart=NULL,alpha=NULL,lamda=NULL,permute=NULL,sig=NULL,calgo=c("rwr","nbi","netcombo","all")){

    A <- t(dataDT)
    S1 = as.matrix(dataD)
    S2 = as.matrix(dataT)

    ## Create a list where to store the matrices
    perm = list()
    print(permute)
    print(sig)
    ## Set a random seed
    set.seed(12345)

    if (calgo=="nbi"){
        P1 <- nbiNet(A,alpha=alpha, lamda=lamda,  s1=S1, s2=S2,format = "matrix")
        for ( i in 1:permute){

            ## randomize the matrices
            B <- A[sample(nrow(A)),sample(ncol(A))]
            S1 <- S1[sample(nrow(S1)),sample(ncol(S1))]
            S2 <- S2[sample(nrow(S2)),sample(ncol(S2))]

            R1 <- nbiNet(B,alpha=alpha, lamda=lamda,  s1=S1, s2=S2,format = "matrix")
            perm[[i]] <- R1
        }

        mean_mat <- apply(simplify2array(perm), 1:2, mean)
        sd_mat  <-  apply(simplify2array(perm), 1:2, sd)

        ## Comput the Z-score of matrix
        Z <- (P1 - mean_mat)/sd_mat
        Z[is.nan(Z)] = 0
        ## Compute the significance score
        sigNetwork <- 2*(pnorm(-abs(Z)))
        results <- getPermScores(sigNetwork,A,sig)
        colnames(results)[1] <- "Drugs"
        colnames(results)[2] <- "Target"
        colnames(results)[3] <- "p-value"
        return(results)
    } else if (calgo=="rwr"){

        g1 = graph.incidence(A)
        Q = biNetwalk(g1,s1=S1,s2=S2,normalise="laplace",dataSeed=NULL,restart=restart,verbose=T)

        for ( i in 1:permute){

            ## randomize the matrices
            B <- A[sample(nrow(A)),sample(ncol(A))]
            S1 <- S1[sample(nrow(S1)),sample(ncol(S1))]
            S2 <- S2[sample(nrow(S2)),sample(ncol(S2))]
            g2 = graph.incidence(B)
            R1 <- biNetwalk(g2,s1=S1,s2=S2,normalise="laplace",dataSeed=NULL,restart=restart,verbose=T)
            perm[[i]] <- R1
        }

        mean_mat <- apply(simplify2array(perm), 1:2, mean)
        sd_mat  <-  apply(simplify2array(perm), 1:2, sd)

        ## Comput the Z-score of matrix
        Z <- (Q - mean_mat)/sd_mat
        Z[is.nan(Z)] = 0
        ## Compute the significance score
        sigNetwork <- 2*(pnorm(-abs(Z)))
        results <- getPermScores(sigNetwork,A,sig)
        colnames(results)[1] <- "Drugs"
        colnames(results)[2] <- "Target"
        colnames(results)[3] <- "p-value"
        return(results)
    }


}

## Get the GO data (BP, CC , MF) for the given set of genes
getGOdata <- function(geneID,level) {
    gbp <- groupGO(gene     = as.character(geneID),
                   OrgDb = 'org.Hs.eg.db',
                   ont      = "BP",
                   level    = as.integer(level),
                   readable = TRUE)
    data <- summary(gbp)
    bp <- data[data$Count > 0,]
    bp$type <- "BP"

    gcc <- groupGO(gene     = as.character(geneID),
                   OrgDb = 'org.Hs.eg.db',
                   ont      = "CC",
                   level    = as.integer(level),
                   readable = TRUE)
    data <- summary(gcc)
    cc <- data[data$Count > 0,]
    cc$type <- "CC"
    gmf <- groupGO(gene     = as.character(geneID),
                   OrgDb = 'org.Hs.eg.db',
                   ont      = "MF",
                   level    = as.integer(level),
                   readable = TRUE)
    data <- summary(gmf)
    mf <- data[data$Count > 0,]
    mf$type <- "MF"
    fgo <- rbind(bp,cc,mf)
    return (fgo)
}

get.ppi <- function(proteins,conf){

    dg <- info.ppi(proteins,conf)


    # netresult <- visNetwork(dg$nodeData, dg$edgeList,background = "lightblue",maxVelocity = 10,minVelocity = 1.0) %>% visIgraphLayout(type= "full",randomSeed = 123,layout = "layout.kamada.kawai") %>% visOptions(selectedBy = "group",highlightNearest = list(enabled = T, hover = T), nodesIdSelection = T) %>%
    #             visLayout(improvedLayout = TRUE) %>% visInteraction(navigationButtons = TRUE,tooltipDelay = 0) %>% visNodes(size = 25,scaling=list(min=20),borderWidth = 2,font=list(size=18),color = list(border = "black",highlight = "yellow")) %>% visLegend(width=0.1) %>% visLayout(randomSeed = 1999,improvedLayout = TRUE) %>% visEvents(click = "function(nodeData) {Shiny.onInputChange('gene_nodeID', nodeData);
    #                                                                                                                                                                                                                                                                                                                                      ;}")
    netresult <- visNetwork(dg$nodeData, dg$edgeList,background = "lightblue",maxVelocity = 10,minVelocity = 1.0)  %>% visOptions(highlightNearest = list(enabled = T, algorithm="hierarchical",hover = T), nodesIdSelection = T) %>%
        visGroups(groupname = "Drug", shape = "icon", icon = list(code = "f055", size = 225, color = "red")) %>% visPhysics(solver = "forceAtlas2Based",forceAtlas2Based = list(gravitationalConstant = -1550,centralGravity=0.003,avoidOverlap=0.88,springLength=50,springConstant = 0.002),stabilization = list(iterations = 200, enabled = TRUE)) %>% visNodes(size = 150,scaling=list(min=20),borderWidth = 2 ,font=list(size=175),color = list(border = "black",highlight = "yellow")) %>% visLayout(improvedLayout = TRUE) %>% visInteraction(navigationButtons = TRUE,tooltipDelay = 0) %>% visEdges(smooth=FALSE) %>% visLegend(width=0.1) %>% visLayout(randomSeed = 1999,improvedLayout = TRUE) %>% visEvents(click = "function(nodeData) {Shiny.onInputChange('gene_nodeID', nodeData);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ;}") %>% addFontAwesome() %>% visExport()

    invisible(netresult)
}

info.ppi <- function(proteins,conf,dataid){

    #myset <- data.frame()


    f <- construct.subgraph(data=as.data.frame(interactome_data),input=as.data.frame(proteins),w=conf,hierarchy=1)

    #f <- construct_subgraph(data=data,input,w=0.99,hierarchy=1)
    f$value <- 1
    fset <- as.data.frame(f)

    #pnodes <- unique(c(as.character(fset$from), as.character(fset$to)))
    #nodelist <- paste("'",as.character(pnodes),"'",collapse =", ",sep ="")

    # if(is.null(dataid)){
    #     return(NULL)
    # }

        # if(dataid == 'show drugs'){
        #
        #     ###############################################
        #     ## Get the drug names associated with Proteins
        #     ##############################################
        #     drugSet <- get.drugs(pnodes)
        #     fset$weight <- NULL
        #     fset <- rbind(fset,drugSet)
        #
        #     #Q1 <- sprintf("SELECT Approved_Symbol FROM hgnc_data where Approved_Symbol IN",nodelist)
        #     #top <- dbGetQuery(con2,Q1)
        #
        #     nodes <- unique(c(as.character(fset$from), as.character(fset$to)))
        #
        #     dcount <- length(unique(drugSet$to))
        #     totalgroup <- c(top$LOCATION,rep("Drug",dcount))
        #
        #     id=seq(0,length(nodes)-1,1)
        #     label= nodes
        #     nodeData <- data.frame(id ,label,shadow=TRUE,group=totalgroup,stringsAsFactors=FALSE)
        #     print(head(nodeData))
        #
        #     #colnames(fullset)[1] <- "from"
        #     #colnames(fullset)[2] <- "to"
        #     set1  <- fset[(fset$from %in% label) & (fset$to %in% label),]
        #     edgeList <- set1
        #     edgeList$from <- with(nodeData, id[match(edgeList$from, label)])
        #     edgeList$to <- with(nodeData, id[match(edgeList$to,label)])
        #
        #     return(list(nodeData = nodeData , edgeList=edgeList))
        #
        # } else
        #if(dataid == 'show ppi'){

            ##################################################
            ### Show the ppi based
            #################################################
            #Q1 <- sprintf("SELECT Approved_Symbol FROM hgnc_data where Approved_Symbol IN (%s)",nodelist)
            #top <- dbGetQuery(con,Q1)

            nodes <- unique(c(as.character(fset$from), as.character(fset$to)))

            #dcount <- length(unique(drugSet$to))
            #totalgroup <- c(top$LOCATION,rep("Drug",dcount))

            id <- seq(0,length(nodes)-1,1)
            label <- nodes
            #nodeData <- data.frame(id ,label,shadow=TRUE,group=top$LOCATION,stringsAsFactors=FALSE)
            nodeData <- data.frame(id ,label,shadow=TRUE,stringsAsFactors=FALSE)

            set1  <- fset[(fset$from %in% label) & (fset$to %in% label),]
            edgeList <- set1
            edgeList$from <- with(nodeData, id[match(edgeList$from, label)])
            edgeList$to <- with(nodeData, id[match(edgeList$to,label)])

            return(list(nodeData = nodeData , edgeList=edgeList))

       # }



    }

get.drugs <- function(proteins){

    myset <- data.frame()

    # proteins <- pnodes

    for ( p in proteins) {

        Q1 <- sprintf("SELECT GENE_NAME,DRUGBANK_ID FROM ABT_DDP.DRUGBANK_DRUG_TARGET WHERE SPECIES = 'Human' AND GENE_NAME = q'[%s]'",p)
        top <- dbGetQuery(con2,Q1)
        if(dim(top)[1] > 0 ){
            mytop <- unique(cSplit(top, "DRUGBANK_ID", ";", direction = "long"))
            myset <- rbind(myset,mytop)
        }

    }

    myset$value <- 1
    colnames(myset) <- c("from","to","value")
    return(myset)
    rm("myset")

}
get.path <- function(from,to,kpath,pathway=FALSE, pval=0.01, weighted){

    #net<-graph.data.frame(interactome_data[,c("SOURCE","TARGET","weight")],directed=FALSE)

    if(is.null(pval)){

        pval <- 0.01
    }

    if(pathway == "FALSE"){

        dg <- info.paths(from,to,kpath,pathway=FALSE,pval=0,weighted)

        netresult <- visNetwork(dg$nodeData, dg$edgeList,background = "lightblue",maxVelocity = 10,minVelocity = 1.0) %>% visNodes(size = 25,scaling=list(min=20,max=50),borderWidth = 2,font=list(size=18),color = list(border = "black",highlight = "yellow")) %>% visEdges(smooth=FALSE) %>% visLegend(width=0.1) %>%
            visOptions(selectedBy = "group",highlightNearest = list(enabled = TRUE ,algorithm="hierarchical",hover=TRUE) ,nodesIdSelection = TRUE)%>% visPhysics(solver = "barnesHut", barnesHut = list(gravitationalConstant = -1500,springLength = 100, springConstant = 0.005),stabilization = list(iterations = 100, enabled = TRUE)) %>% visLayout(improvedLayout = TRUE) %>% visInteraction(navigationButtons = TRUE,tooltipDelay = 0) %>% visEvents(click = "function(nodeData) {Shiny.onInputChange('path_nodeID', nodeData);
                                                                                                                                                                                                                                                                                                                                                                                                                                                           ;}") %>% visExport()
        return(list(netresult = netresult))

} else {

    ## Perform pathway Enrichment along with PPI

    dg <- info.paths(from,to,kpath,pathway=TRUE,pval=pval,weighted)

    netresult <- visNetwork(dg$nodeData, dg$edgeList,background = "lightblue",maxVelocity = 10,minVelocity = 1.0)  %>% visOptions(highlightNearest = list(enabled = T, algorithm="hierarchical",hover = T), nodesIdSelection = T) %>%
        visPhysics(solver = "forceAtlas2Based",forceAtlas2Based = list(gravitationalConstant = -1600,centralGravity=0.001,avoidOverlap=0.9,springLength=50,springConstant = 0.002),stabilization = list(iterations = 100, enabled = TRUE)) %>% visNodes(size = 250,scaling=list(min=120,max=300),borderWidth = 2 ,font=list(size=150),color = list(border = "black",highlight = "yellow")) %>% visLayout(improvedLayout = TRUE) %>% visInteraction(navigationButtons = TRUE,tooltipDelay = 0) %>% visEdges(smooth=FALSE) %>% visLegend(width=0.1) %>% visLayout(randomSeed = 1999,improvedLayout = TRUE) %>% visExport()
    #return(netresult)
    return(list(netresult = netresult , pathway = dg$path))
}

    rm("netresult")
}

get.sub <- function(genes,dataset="ConsensusPath") {
  
  print (genes)
  if (dataset == "ConsensusPath"){
    
    graph <- Weight_ConsensusPath
    
    new_graph <- induced.subgraph(graph,which(V(graph)$name %in% genes))
    
  } else if (dataset == "String"){
    
    
    graph <- Weight_string
    
    new_graph <- induced.subgraph(graph,which(V(graph)$name %in% genes))
    
  }
  
  subNet <- unique(get.data.frame(new_graph))
  
  
  subNet$value <- 1
  fset <- as.data.frame(subNet)
  
  #pnodes <- unique(c(as.character(fset$from), as.character(fset$to)))
  #nodelist <- paste("'",as.character(pnodes),"'",collapse=", ",sep="")
  
  
  nodes <- unique(c(as.character(fset$from), as.character(fset$to)))

  
  id=seq(0,length(nodes)-1,1)
  label= nodes
  nodeData <- data.frame(id ,label,shadow=TRUE,stringsAsFactors=FALSE)
  set1  <- fset[(fset$from %in% label) & (fset$to %in% label),]
  edgeList <- set1
  edgeList$from <- with(nodeData, id[match(edgeList$from, label)])
  edgeList$to <- with(nodeData, id[match(edgeList$to,label)])
  
  netresult <- visNetwork(nodeData, edgeList,background = "lightblue",maxVelocity = 10,minVelocity = 1.0)  %>% visOptions(selectedBy = "group",highlightNearest = list(enabled = T, algorithm="hierarchical",hover = T), nodesIdSelection = T) %>% 
    visPhysics(solver = "forceAtlas2Based",forceAtlas2Based = list(gravitationalConstant = -1550,centralGravity=0.003,avoidOverlap=0.88,springLength=50,springConstant = 0.002),stabilization = list(iterations = 200, enabled = TRUE)) %>% visNodes(size = 80,scaling=list(min=100,max=175),borderWidth = 2 ,font=list(size=175),color = list(border = "black",highlight = "yellow")) %>% visLayout(improvedLayout = TRUE) %>% visInteraction(navigationButtons = TRUE,tooltipDelay = 0) %>% visEdges(smooth=FALSE) %>% visLegend(width=0.1) %>% visLayout(randomSeed = 1999,improvedLayout = TRUE) %>% visEvents(click = "function(nodeData) {Shiny.onInputChange('gene_nodeID', nodeData);;}") %>% addFontAwesome() %>% visExport()
  
  return(netresult)
  
  
  rm("netresult")
  
}
## Getting Reactome Pathways and PPI information

info.paths <- function(from,to,kpath,pathway=FALSE,pval=0.01,weighted='Weighted PPI'){

    if (weighted == 'Weighted PPI'){

        net <- graph.data.frame(interactome_data[,c("SOURCE","TARGET","weight")],directed=FALSE)
    } else if ((weighted == 'Unweighted PPI')){
        interactome_data$weight <- 1
        net <- graph.data.frame(interactome_data[,c("SOURCE","TARGET")],directed=FALSE)

    }
    print (weighted)
    #print (tail(net))

    k_path <- k.shortest.paths(net, from = from, to = to, k = kpath)
    elist <- data.frame()

    for ( i in 1:length(k_path)){
        vseq <- k_path[[i]]$vert[[1]]$name
        cvert <- length(vseq) - 1
        for ( j in 1:cvert) {
            edge <- cbind(vseq[j],vseq[j+1])
            elist <- rbind(edge,elist)
        }
    }

    elist <- unique(elist)
    elist$value <- 1
    colnames(elist) <- c("from","to","value")

    ## Search a Direct interaction if exists and connect to the paths. So a SQL query.

    #Q2 <- sprintf("SELECT HGNC_SOURCE,HGNC_TARGET,CONFIDENCE from ABT_DDP.ABBV_GENE_PPI where HGNC_SOURCE = '%s' AND HGNC_TARGET = '%s'",from,to)

    #top2 <- dbGetQuery(con2,Q2)
    top2 <- filter(interactome_data, SOURCE == from & TARGET == to)
    colnames(top2) <- c("from","to","value")

    #Q3 <- sprintf("SELECT HGNC_SOURCE,HGNC_TARGET,CONFIDENCE from ABT_DDP.ABBV_GENE_PPI where HGNC_SOURCE = '%s' AND HGNC_TARGET = '%s'",to,from)

    #top3 <- dbGetQuery(con2,Q3)


    top3 <- filter(interactome_data, TARGET == from & SOURCE == to)

    colnames(top3) <- c("from","to","value")


    elist <- rbind(elist,top2,top3)
    #elist <- rbind(elist,top3)

    if (pathway == FALSE){

        nodes <- unique(c(as.character(elist$from), as.character(elist$to)))
        #nodelist <- paste("'",as.character(nodes),"'",collapse=", ",sep="")


        ## Get Location of the proteins
        #Q1 <- sprintf("SELECT HGNC,LOCATION from ABT_DDP.ABBV_GENE_META where HGNC IN (%s)",nodelist)

        #top1 <- dbGetQuery(con2,Q1)

        id=seq(0,length(nodes)-1,1)
        label=nodes
        #nodeData <- data.frame(id ,label,shadow=TRUE,group=top1$LOCATION,stringsAsFactors=FALSE)
        nodeData <- data.frame(id ,label,shadow=TRUE,stringsAsFactors=FALSE)
        nodeData$value <- 15

        nodeData$value[nodeData$label==from] <- 40
        nodeData$value[nodeData$label==to] <- 40

        edgeList <- elist
        edgeList$from <- with(nodeData, id[match(edgeList$from, nodes)])
        edgeList$to <- with(nodeData, id[match(edgeList$to,nodes)])




    } else if(pathway == TRUE){


        pnodes <- unique(c(as.character(elist$from), as.character(elist$to)))
        #pnodelist <- paste("'",as.character(pnodes),"'",collapse=", ",sep="")

        ## Get location of the proteins
        #Q1 <- sprintf("SELECT HGNC,LOCATION from ABT_DDP.ABBV_GENE_META where HGNC IN (%s)",pnodelist)
        #top1 <- dbGetQuery(con2,Q1)


        mapTab <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = pnodes, mart = human, uniqueRows=FALSE)
        geneID <- unique(mapTab$entrezgene)
        x <- enrichPathway(gene = as.character(geneID),organism = "human",pvalueCutoff = 0.1,pAdjustMethod = "BH", readable = T)
        path <- summary(x)
        p <- path[path$Count > 0,]
        p <- p[p$p.adjust <= pval,]
        pathVars <- p[,c('Description','geneID')]
        colnames(pathVars) <- c("Pathway","gene")
        print(head(pathVars))

        piv <- pathVars %>% mutate(gene = strsplit(as.character(gene),"/")) %>% unnest(gene)
        piv$value <- 1
        colnames(piv) <- c("to","from","value")

        elist <- rbind(elist,piv)
        #totalPath <- length(unique(piv$to))
        #group <- c(top1$LOCATION,rep("Pathway",totalPath))

        nodes <- unique(c(as.character(elist$from), as.character(elist$to)))
        id=seq(0,length(nodes)-1,1)
        label=nodes
        nodeData <- data.frame(id ,label,shadow=TRUE,stringsAsFactors=FALSE)
        nodeData$value <- 30

        nodeData$value[nodeData$label==from] <- 50
        nodeData$value[nodeData$label==to] <- 50

        edgeList <- elist
        edgeList$from <- with(nodeData, id[match(edgeList$from, nodes)])
        edgeList$to <- with(nodeData, id[match(edgeList$to,nodes)])

        # netresult
        #return(list(nodeData = nodeData , edgeList=edgeList))


    }
    return(list(nodeData = nodeData , edgeList=edgeList, pathway = path))

}

## Add attributes to the vertex.
add.vertex.attri <- function(graph=NULL,input)
{
    if(!is.igraph(graph)){
        stop("Not an igraph object")
    }
    if(!is.data.frame(input)){
        input <- as.data.frame(input)
    }
    coln <- ncol(input)
    if(coln > 1) {
        colname <- colnames(input)
        index <- match(V(graph)$name,input[,1],nomatch=0)
        for(i in 2:coln) {
            graph <- set.vertex.attribute(graph,colname[i],index=as.character(input[index,1]),
                                        value=input[index,i])
        }
    }
    return(graph)
}


## PPI based subgraph construction code

construct.subgraph <- function(input,data,w=1,hierarchy=1){


    data <- data[data$weight >= w,]
    if(!missing(input)){
        if( !is.data.frame(input)){
            stop("Not a data frame")
        }
    }

    ##import the data

    net<-graph.data.frame(data[,c("SOURCE","TARGET","weight")],directed=FALSE)

    if(missing(input)){
        return(net)
    }

    ##  match
    index <- match(input[,1],V(net)$name)
    index <- index[!is.na(index)]
    ##  create a  sub network
    graph <- induced.subgraph(graph=net,unlist(neighborhood(graph=net,order=hierarchy,nodes=index)))
    ##  vertex.hierarchy
    V(graph)$vertex.hierarchy <- rep(hierarchy,vcount(graph))
    if(hierarchy >= 1){
        for(i in (hierarchy-1):0){
            index <- unlist(neighborhood(graph=graph,order=i,nodes=which(V(graph)$name %in% input[,1])))
            V(graph)$vertex.hierarchy[index] <- i
        }
    }
    ##  add vertex attributes
    graph <- add.vertex.attri(graph=graph,input=input)

    gframe <- unique(get.data.frame(graph))
    return(gframe)
}

# find k shortest paths
k.shortest.paths <- function(graph, from, to, k){
    # first shortest path
    k0 <- get.shortest.paths(graph,from,to, output='both')

    # number of currently found shortest paths
    kk <- 1

    # list of alternatives
    variants <- list()

    # shortest variants

    # shortest variants
    shortest.variants <- list(list(g=graph, path=k0$epath, vert=k0$vpath, dist=shortest.paths(graph,from,to)))

    # until k shortest paths are found
    while(kk < k){
        # take last found shortest path
        last.variant <- shortest.variants[[length(shortest.variants)]]

        # calculate all alternatives
        variants <- calculate.variants(variants, last.variant, from, to)

        # find shortest alternative
        sp <- select.shortest.path(variants)

        # add to list, increase kk, remove shortest path from list of alternatives
        shortest.variants[[length(shortest.variants)+1]] <- list(g=variants[[sp]]$g, path=variants[[sp]]$variants$path, vert=variants[[sp]]$variants$vert, dist=variants[[sp]]$variants$dist)
        kk <- kk+1
        variants <- variants[-sp]
    }

    return(shortest.variants)
}

# found all alternative routes
calculate.variants <- function(variants, variant, from, to){
    # take graph from current path
    g <- variant$g

    # iterate through edges, removing one each iterations
    for (j in unlist(variant$path)){
        newgraph <- delete.edges(g, j) # remove adge
        sp <- get.shortest.paths(newgraph,from,to, output='both') # calculate shortest path
        spd <- shortest.paths(newgraph,from,to) # calculate length
        if (spd != Inf){ # the the path is found
            if (!contains.path(variants, sp$vpath)) # add to list, unless it already contains the same path
            {
                variants[[length(variants)+1]] <- list(g=newgraph, variants=list(path=sp$epath, vert=sp$vpath, dist=spd))
            }
        }
    }

    return(variants)
}

# does a list contain this path?
contains.path <- function(variants, variant){
    return( any( unlist( lapply( variants, function(x){ identical(unlist(x$variant$vert),unlist(variant)) } ) ) ) )
}

# which path from the list is the shortest?
select.shortest.path <- function(variants){
    return( which.min( unlist( lapply( variants, function(x){x$variants$dist} ) ) ) )
}
