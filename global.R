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
# list.of.packages <- c("devtools","DBI","RSQLite","doParallel","foreach","shinyGridster","shinyBS","netpredictor","lpbrim","ShinySky","shinyjs","shinythemes","reshape2","rlist","htmltools","igraph","gdata","shiny","data.table","visNetwork","DBI","RSQLite") # replace xx and yy with package names
# 
# packages.auto <- function(x) { 
#     x <- as.character(substitute(x)) 
#     if(isTRUE(x %in% .packages(all.available=TRUE))) { 
#         eval(parse(text = sprintf("require(\"%s\")", x)))
#     } else { 
#         #update.packages(ask= FALSE) #update installed packages.
#         eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
#     }
#     if(isTRUE(x %in% .packages(all.available=TRUE))) { 
#         eval(parse(text = sprintf("require(\"%s\")", x)))
#     } else {
#         
#         devtools::install_github("abhik1368/netpredictor")
#         devtools::install_github("AnalytixWare/ShinySky")
#         devtools::install_github("PoisotLab/lpbrim")
#         devtools::install_github("shiny-gridster", "wch")
#         devtools::install_github("ebailey78/shinyBS")
# 
#         
#     }
# }
# 
# packages.auto(list.of.packages)
#################################################################################

#################################################################################
## Load the libraries

library(netpredictor)
library(igraph)
library(shiny)
require(shinysky)
library(reshape2)
library(lpbrim)
library(rlist)
set.seed(12345)
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
