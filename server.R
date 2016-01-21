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
## Load the libraries

library(netpredictor)
library(igraph)
library(gdata)
library(shiny)
require(shinysky)
#library(rcytoscapejs)
library(data.table)
## Server for netpredictor
source('global.R')
library(networkD3)
library(visNetwork)
<<<<<<< HEAD
library(DBI)
library(RSQLite)
set.seed(12345)

=======
#################################################################################
>>>>>>> origin/master

shinyServer( function(input, output,session) {
    
    observe({
        if(input$gofind){
            updateTabsetPanel(session, "bar", selected = "Start Prediction")
        }
    })
    
    observe({
        if(input$start){
            updateTabsetPanel(session, "datatabs", selected = "Prediction Results")
        }
    })
    
    
## Get the properties         
 prop <- reactive({
     
     if (input$netproperty <= 0){
         return(NULL)
     } 
     result <- isolate({ 
         input$netproperty
         tryCatch ({
         
         if(input$data_input_type=="example"){
             
             if(input$datasets == "Enzyme"){
                 load("Enzyme.rda")
                 data <- t(adjm)
                 props <- getProp(data)
                 props
             } else if (input$datasets == "GPCR"){
                 load("GPCR.rda")
                 data <- t(adjm)
                 props <- getProp(data)
                 props
             } else if (input$datasets == "Nuclear Receptor"){
                load("Nuclear Receptor.rda")
                 data <- t(adjm)
                 props <- getProp(data)
                 props
             } else if (input$datasets == "Ion Channel") {
                load("Ion Channel.rda")
                 data <- t(adjm)
                 props <- getProp(data)
                 props
             } 
        } else if (input$data_input_type=="custom"){
            if (is.null(input$dt_file))
                return(NULL)
            inFile <- input$dt_file
            dataDT <- as.matrix(read.csv(inFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
            props <- getProp(dataDT)
            print (props)
            props
         }},
         error = function(err) {
             print(paste("Select Example or custom data set"))
         })
     })

    })
 
## Get the modules of the bipartite network
totModules <- reactive({
     
    if (input$mods <= 0){
             return(NULL)
    } 
    
     input$mods
     results <- isolate({
         input$mods
         if(input$data_input_type=="example"){
             dset <- input$datasets
             exdata <- paste(dset,".rda",sep="")
             load(exdata)
             dm <-as.matrix(adjm)
             mod <- getMod(dm)
             mod
         } else if(input$data_input_type=="custom") {
             if (is.null(input$dt_file))
                 return(NULL)
             
             DTFile <- input$dt_file
             dataDT <- as.matrix(read.csv(DTFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
             mod <- getMod(dataDT)
             mod
         }
         
     })
     results
 })

## This for the drugs 
topDrugs <- reactive({
     
    if (input$netproperty <= 0){
        return(NULL)
    } 
        # Use isolate() to avoid dependency on input$obs
        results <- isolate({
            input$netproperty
            if(input$data_input_type=="example"){
                dset <- input$datasets
                exdata <- paste(dset,".rda",sep="")
                load(exdata)
                dt <- t(adjm)
                topRows <- data.frame(rowSums(dt))
                colnames(topRows)[1] <- "count"
                topRows$Drugs <- row.names(topRows)
                res <- topRows[order(-topRows$count),][1:15,]
                rownames(res)<- NULL
                res 
            } else if(input$data_input_type=="custom") {
                if (is.null(input$dt_file))
                    return(NULL)
                
                DTFile <- input$dt_file
                dataDT <- as.matrix(read.csv(DTFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
                topRows <- data.frame(colSums(dataDT))
                colnames(topRows)[1] <- "count"
                topRows$Drugs <- row.names(topRows)
                res <- topRows[order(-topRows$count),][1:15,]
                rownames(res)<- NULL
                res
            }
            
            })
        results
 })
 
topProteins <- reactive({
    
    if (input$netproperty <= 0){
        return(NULL)
    } 
    
    # Use isolate() to avoid dependency on input$obs
    results <- isolate({
        input$netproperty
        if(input$data_input_type=="example"){
            dset <- input$datasets
            exdata <- paste(dset,".rda",sep="")
            load(exdata)
            dt <- t(adjm)
            ## This for the drugs
            topRows <- data.frame(colSums(dt))
            colnames(topRows)[1] <- "count"
            topRows$Proteins <- row.names(topRows)
            res <- topRows[order(-topRows$count),][1:15,]
            rownames(res)<- NULL
            res 
        } else if(input$data_input_type=="custom") {
            if (is.null(input$dt_file))
                return(NULL)
            
            DTFile <- input$dt_file
            dataDT <- as.matrix(read.csv(DTFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
            topRows <- data.frame(rowSums(dataDT))
            colnames(topRows)[1] <- "count"
            topRows$Proteins <- row.names(topRows)
            res <- topRows[order(-topRows$count),][1:15,]
            rownames(res)<- NULL
            res
        }
        
    })
    results
})




## Output properties table 
output$prop_table <- renderTable({
    prop()
     
 })
 
## Output count distribution of drugs
output$countDrugs <- renderTable({
    topDrugs()
})

## Output count distribution of proteins
output$countProteins <- renderTable({
    topProteins()
})

## Get the count of modules

tot_modules <- reactive({
    mod <- totModules()
    len <- length(mod)
    modules <- c()
    for (i in 1:len){
        
        modules[i] <- paste("module",toString(i),sep="")
    }
   return(modules)
 })

## Get the list of modules
full_modules <- reactive({
    mod <- totModules()
    return(mod)
})

## Update dynamically the list of modules
output$modules = renderUI({
    
    selectInput('module', 'Modules', tot_modules())
})


# Output the data table of modules
data_table <- reactive({
    # If missing input, return to avoid error later in function
    if(is.null(input$module))
        return()
    
    # Get the module
    dat <- input$module
    
    if (dat %in%  tot_modules())
        modIndx <- match(dat, tot_modules())
    
    modMat <- full_modules()
    mod <- modMat[[modIndx]]
    gMod <- graph.incidence(mod)
    el <- get.edgelist(gMod)
    colnames(el) <- c("column1","column2")
    return (el)
    
})

## Out the module list on the data table
output$data_table <- renderDataTable( {
    data_table()
},options = list(pageLength = 10,
    autoWidth = TRUE,
    columnDefs = list(list(width = '50px', targets = "_all")))
)


modnetwork <- reactive({
    
    if (input$shownet <= 0){
        return(NULL)
    }
    
    netResult <- isolate({
        input$shownet
        if(is.null(input$module))
            return()
        
        # Get the module
        dat <- input$module
        
        if (dat %in%  tot_modules())
            modIndx <- match(dat, tot_modules())
        
        modMat <- full_modules()
        mod <- modMat[[modIndx]]
        gMod <- graph.incidence(mod)
        el <- data.frame(get.edgelist(gMod))
        colnames(el) <- c("column1","column2")
        nodes <- unique(c(as.character(el$column1), as.character(el$column2)))
        group1 <- length(unique(el$column1))
        group2 <- length(unique(el$column2))
        id=seq(0,length(nodes)-1,1)
        label=nodes
        print (dim(el))
        nodeData <- data.frame(id ,label,group=c(rep("Drugs",group1),rep("Proteins",group2)),stringsAsFactors=FALSE)
        colnames(el)[1] <- "from"
        colnames(el)[2] <- "to"
        edgeList <- el[, c("from","to")]
        edgeList$from <- with(nodeData, id[match(edgeList$from, nodes)])
        edgeList$to <- with(nodeData, id[match(edgeList$to,nodes)])
        #edgeList$dashes <- ifelse(mynet$type == "True Interactions",FALSE,TRUE)
        netresult <- visNetwork(nodeData, edgeList,width = "70%") %>% visNodes(size = 25) %>% visInteraction(navigationButtons = TRUE,tooltipDelay = 0) %>% visLegend() %>% visNodes(scaling=list(min=20),font=list(size=24)) %>%
            visOptions(selectedBy = "group", nodesIdSelection = TRUE,highlightNearest = TRUE)%>% visPhysics(solver = "barnesHut",barnesHut = list(gravitationalConstant = -1000,avoidOverlap=0.5,springLength=250))
        
        
    })
    netResult
})


output$moduleplot <- renderVisNetwork({
    modnetwork()
})

 observe({
     if (input$start == 0) 
         return()
     showshinyalert(session, "shinyalert1", paste("Computation in progress", "success"), 
                    styleclass = "success")
 })

 
 Result <- reactive({
     
     
     if (input$start <= 0){
         return(NULL)
     } 
     
     # Take a dependency on input$start
     input$start
     
     if(input$data_input_type=="example"){
         # Use isolate() to avoid dependency on input$obs
         results <- isolate({
    
         dset <- input$datasets
         exdata <- paste(dset,".rda",sep="")
         
         if (input$algorithm_typeI == "heat"){
             heat_alpha = input$heatAlpha
             heat_lamda = input$heatLambda
             results <- getheat(exdata,heat_alpha,heat_lamda)
             results
         } else if (input$algorithm_typeI == "nbi"){
             nbi_alpha = input$nbiAlpha
             nbi_lamda = input$nbiLambda
             results <- getnbi(exdata,nbi_alpha,nbi_lamda)
             results
         } else if (input$algorithm_typeI == "rwr"){
             restart = input$rwr_restart
             results <- getrwr(exdata,restart)
             results
         } else if (input$algorithm_typeI == "netcombo"){
             restart = input$nc_restart
             nc_alpha = input$ncAlpha
             nc_lamda = input$ncLambda
             results <- getcombo(exdata,nc_alpha,nc_lamda,restart)
             results
         } 
     })
     
      results
     } else if (input$data_input_type=="custom"){
         if (is.null(input$dt_file))
             return(NULL)
         if (is.null(input$drug_file))
             return(NULL)
         if (is.null(input$target_file))
             return(NULL)
         
         results <- isolate({ 
             
             DTFile <- input$dt_file
             dataDT <- as.matrix(read.csv(DTFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
             
             DFile <- input$drug_file
             dataD <- as.matrix(read.csv(DFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
             
             TFile <- input$target_file
             dataT <- as.matrix(read.csv(TFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
             
             if (input$algorithm_typeII == "heat"){
                 heat_alpha = input$cheatAlpha
                 heat_lamda = input$cheatLambda
                 results <- getCustomHeat(dataDT,heat_alpha,heat_lamda)
                 results
             } else if (input$algorithm_typeII == "nbi"){
                 nbi_alpha = input$cnbiAlpha
                 nbi_lamda = input$cnbiLambda
                 results <- getCustomNBI(dataDT,dataD,dataT,nbi_alpha,nbi_lamda)
                 results
             } else if (input$algorithm_typeII == "rwr"){
                 restart = input$crwr_restart
                 results <- getCustomRWR(dataDT,dataD,dataT,restart)
                 results
             } else if (input$algorithm_typeII == "netcombo"){
                 restart = input$cnc_restart
                 nc_alpha = input$cncAlpha
                 nc_lamda = input$cncLambda
                 results <- getCustomCombo(dataDT,dataD,dataT,nc_alpha,nc_lamda,restart)
                 results
             } 
         })
         results
     }
     
     
     
 })
 
 ## Get the output Result to data table output
 
 output$Result <- renderDataTable(Result(),options = list(pageLength = 10))
 
 
 output$downloadResult <- downloadHandler(
     
     filename = function() { paste("Results.txt") },
     
     content = function(file) {
         write.table(Result(), file, row.names=FALSE, quote=FALSE, sep="\t")
         
     })
 

 
#  output$networkplot <- renderForceNetwork({
#      
#      if (input$start <= 0){
#          return(NULL)
#      } 
#      input$start
#      netResult <- isolate({
#          mynet <- Result()
#          nodes <- unique(c(mynet$Proteins, mynet$Drugs))
#          group1 <- length(unique(mynet$Proteins))
#          group2 <- length(unique(mynet$Drugs))
#          id = nodes
#          name = nodes
#          #nodeData <- data.frame(name,stringsAsFactors=FALSE)
#          nodeData <- data.frame(name,group=c(rep(10,group1),rep(25,group2)),stringsAsFactors=FALSE)
#          colnames(mynet)[1] <- "source"
#          colnames(mynet)[2] <- "target"
#          colnames(mynet)[3] <- "value"
#          edgeList <- mynet[, c("source","target","value")]
#          #edgeData <- edgeList
#          
#          nodesID <- data.frame(nodes,id=seq(0,length(nodes)-1,1))
#          edgeList$source <- with(nodesID, id[match(edgeList$source, nodes)])
#          edgeList$target <- with(nodesID, id[match(edgeList$target,nodes)])
#          nodeData$name <- factor(nodeData$name)
#          
#          
#          forceNetwork(Links = edgeList, Nodes = nodeData, Source = "source",
#                       Target = "target", NodeID = "name",zoom = TRUE,Value="value",
#                       Group = "group",fontFamily="verdana",linkColour = "#afafaf",opacity = 0.8,legend=T)
#          #cyNetwork <- createCytoscapeJsNetwork(nodeData, edgeData)
#          #rcytoscapejs(nodeEntries=cyNetwork$nodes, edgeEntries=cyNetwork$edges)
#      })
#     netResult
#  })
#  
 
 output$networkplot <- renderVisNetwork({
     
     if (input$start <= 0){
         return(NULL)
     }
     input$start
     netResult <- isolate({
         mynet <- Result()
         nodes <- unique(c(as.character(mynet$Proteins), as.character(mynet$Drugs)))
         group1 <- length(unique(mynet$Proteins))
         group2 <- length(unique(mynet$Drugs))
         id=seq(0,length(nodes)-1,1)
         label=nodes
         
         nodeData <- data.frame(id ,label,group=c(rep("Proteins",group1),rep("Drugs",group2)),stringsAsFactors=FALSE)
         colnames(mynet)[1] <- "from"
         colnames(mynet)[2] <- "to"
         edgeList <- mynet[, c("from","to")]
         edgeList$from <- with(nodeData, id[match(edgeList$from, nodes)])
         edgeList$to <- with(nodeData, id[match(edgeList$to,nodes)])
         edgeList$dashes <- ifelse(mynet$type == "True Interactions",FALSE,TRUE)
         netresult <- visNetwork(nodeData, edgeList,height="700px",width = "100%") %>% visInteraction(navigationButtons = TRUE,tooltipDelay = 0) %>% visNodes(size=25) %>% visNodes(scaling=list(min=20),font=list(size=24)) %>%
             visOptions(selectedBy = "group", nodesIdSelection = TRUE,highlightNearest = TRUE) %>% visPhysics(stabilization=FALSE) %>% visPhysics(solver = "barnesHut",barnesHut = list(gravitationalConstant = -1000,avoidOverlap=0.5,springLength=250))
         
 
     })
     netResult
 })
 
 
 output$graphResult <- downloadHandler(
     
     filename = function() { paste("graph.gml") },
     
     content = function(file) {
         netResult <- Result()
         g<-graph.data.frame(netResult[,1:2],directed=FALSE)
         
         ## Set the edge values
         g <- set.edge.attribute(g, "weight", value=netResult[,3])
         saveGML(g,file ,"netresult")
         
     })

 
 ## For the Statistical Analysis
 
 advancedResult <- reactive({ 
     
     input$submit
     if(input$data_input_type=="example"){
         # Use isolate() to avoid dependency on input$obs
         adresults <- isolate({
             
             dset <- input$datasets
             exdata <- paste(dset,".rda",sep="")
             #print(exdata)
             load(exdata)
             rl = input$relinks
             nT = input$freqT
             if (input$predMetrics == "nbi"){
                 nbi_alpha = input$pdnbiAlpha
                 nbi_lamda = input$pdnbiLambda
                 results <- data.frame(net.perf(adjm,protSim,compSim,alpha=nbi_alpha,lamda=nbi_lamda,relinks =rl,numT=nT,Calgo="nbi"))
                 results$num_links_remove = rl
                 results$Freq_Association = nT
                 results$method = "NBI"
                 results
             } else if (input$predMetrics == "rwr"){
                 rt = input$pdrwrRestart
                 results <- data.frame(net.perf(adjm,protSim,compSim,restart=rt, relinks =rl,numT=nT,Calgo="rwr"))
                 results$num_links_remove = rl
                 results$Freq_Association = nT
                 results$method = "RWR"
                 print (head(results))
                 results
                 
             } else if (input$predMetrics == "nc"){
                 nc_rt = input$pdcnc_restart
                 nc_alpha = input$pdcncAlpha
                 nc_lamda = input$pdcncLambda
                 results <- data.frame(net.perf(adjm,protSim,compSim,restart=nc_rt,alpha=nc_alpha,lamda=nc_lamda,relinks =rl,numT=nT,Calgo="netcombo"))
                 results$num_links_remove = rl
                 results$Freq_Association = nT
                 results$method = "NetCombo"
                 results
             } 
         })
     
     
     } else if (input$data_input_type=="custom"){
     
         if (is.null(input$dt_file))
             return(NULL)
         if (is.null(input$drug_file))
             return(NULL)
         if (is.null(input$target_file))
             return(NULL)
         
         adresults <- isolate({ 

             DTFile <- input$dt_file
             dataDT <- as.matrix(read.csv(DTFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
             
             DFile <- input$drug_file
             dataD <- as.matrix(read.csv(DFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
             
             TFile <- input$target_file
             dataT <- as.matrix(read.csv(TFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
             rl = input$relinks
             nT = input$freqT
             
            if (input$predMetrics == "nbi"){
                 nbi_alpha = input$pdnbiAlpha
                 nbi_lamda = input$pdnbiLambda
                 results <- data.frame(net.perf(dataDT,dataT,dataD,alpha=nbi_alpha,lamda=nbi_lamda,relinks =rl,numT=nT,Calgo="nbi"))
                 results$num_links_remove = rl
                 results$Freq_Association = nT
                 results$method = "NBI"
                 results
             } else if (input$predMetrics == "rwr"){
                 rt = input$pdrwrRestart
                 results <- data.frame(net.perf(dataDT,dataT,dataD,restart=rt, relinks =rl,numT=nT,Calgo="rwr"))
                 results$num_links_remove = rl
                 results$Freq_Association = nT
                 results$method = "RWR"
                 print(head(results))
                 results
             } else if (input$predMetrics == "nc"){
                 nc_rt = input$pdcnc_restart
                 nc_alpha = input$pdcncAlpha
                 nc_lamda = input$pdcncLambda
                 results <- data.frame(net.perf(dataDT,dataT,dataD,restart=nc_rt,alpha=nc_alpha,lamda=nc_lamda,relinks =rl,numT=nT,Calgo="netcombo"))
                 results$num_links_remove = rl
                 results$Freq_Association = nT
                 results$method = "NetCombo"
                 results
             } 
         })
     }
         
     })
     
     myresult <- reactiveValues()
     myresult$df <- data.frame()
     newEntry <- observe({
         if(input$submit>0){
             newdata <- isolate(advancedResult())
             isolate(myresult$df <- rbind(myresult$df,newdata))
         }
     })
     
     
     
     
    ## This function for the significance results tab 
    sigResult <- reactive({ 
        
        
        if (input$sigSubmit <= 0){
            return(NULL)
        } 
         input$sigSubmit

         if(input$data_input_type=="example"){
             
             adresults <- isolate({
                 permutation <- input$permute
                 significance <- input$sig
                 dset <- input$datasets
                 exdata <- paste(dset,".rda",sep="")
                 #print(exdata)
                 load(exdata)
                 if (input$sigMetrics == "signbi"){
                     nbi_alpha = input$sgnbiAlpha
                     nbi_lamda = input$sgnbiLambda
                     sigResults <- permTest(adjm,compSim,protSim,alpha=nbi_alpha, lamda=nbi_lamda, permute=permutation, sig=significance,calgo="nbi")
                     sigResults
                 } else if (input$sigMetrics == "sigrwr"){
                     rt = input$sgrwrRestart
                     sigResults <- permTest(adjm,compSim,protSim,restart=rt, perm=permutation, sig=significance, calgo="rwr")
                     sigresults
                 }
             })
         
             
         }  else if(input$data_input_type=="custom"){
             
             if (is.null(input$dt_file))
                 return(NULL)
             if (is.null(input$drug_file))
                 return(NULL)
             if (is.null(input$target_file))
                 return(NULL)
             
             adresults <- isolate({
                 DTFile <- input$dt_file
                 dataDT <- as.matrix(read.csv(DTFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
                 
                 DFile <- input$drug_file
                 dataD <- as.matrix(read.csv(DFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
                 
                 TFile <- input$target_file
                 dataT <- as.matrix(read.csv(TFile$datapath, sep=",", header=TRUE, fill=TRUE,row.names = 1,quote = '"'))
                 
                 permutation <- input$permute
                 significance <- input$sig

                 if (input$sigMetrics == "signbi"){
                     nbi_alpha = input$sgnbiAlpha
                     nbi_lamda = input$sgnbiLambda
                     sigResults <- permTest(dataDT,dataD,dataT,alpha=nbi_alpha, lamda=nbi_lamda, permute=permutation, sig=significance,calgo="nbi")
                     sigResults
                 } else if (input$sigMetrics == "sigrwr"){
                     rt = input$sgrwrRestart
                     sigResults <- permTest(dataDT,dataD,dataT,restart=rt, perm=permutation, sig=significance, calgo="rwr")
                     sigResults
                 }
             })
             
             
         } 
           
        })
     
 ## output table for the navbarmenu tabs    
 output$advTable <-  renderDataTable(myresult$df)
 
 output$downloadadvr <- downloadHandler(
     
     filename = function() { paste("advResults.txt") },
     
     content = function(file) {
         write.table(myresult$df, file, row.names=FALSE, quote=FALSE, sep="\t")
         
     })
 

  output$sigTable <-  renderDataTable({
     sigResult()
     })
 
  output$downloadSig <- downloadHandler(
      
      filename = function() { paste("Results.txt") },
      
      content = function(file) {
          write.table(sigResult(), file, row.names=FALSE, quote=FALSE, sep="\t")
          
      })
  
  
  ## Drugbank Tabs
  
  dResult <- reactive({ 
      if (input$dSearch <= 0){
          return(NULL)
      } 
      
      input$dSearch
      dresults <- isolate({
         if(input$search_type=='drugs'){
             input$dSearch
             drug = input$did
             dname <-  gsub(" ", "", drug, fixed = TRUE)
             con = dbConnect(SQLite(), dbname="drugbank_prediction.db")
             #alltables = dbListTables(con)
             dq <-  sprintf("select DRUGBANK_ID,drugbank_predictions.UNIPROTID,pvalue,outcome,NAME,GENE from drugbank_predictions,target_info where drugbank_predictions.UNIPROTID=target_info.UNIPROTID and drugbank_predictions.DRUGBANK_ID=\'%s\'",dname)
             p1 = dbGetQuery(con,dq)
             dbDisconnect(con)
             p1
         } else if(input$search_type=='proteins'){
              input$dSearch
              protein = input$pid
              pname <-  gsub(" ", "", toupper(protein), fixed = TRUE)
              
              con = dbConnect(SQLite(), dbname="drugbank_prediction.db")
              pq <-  sprintf("select DRUGBANK_ID,GENE,pvalue,outcome from drugbank_predictions, target_info where drugbank_predictions.UNIPROTID = target_info.UNIPROTID and target_info.GENE =\'%s\'",pname)
              t1 = dbGetQuery(con,pq)
              dbDisconnect(con)
              t1
          }
  
      })
        #dresults
  })

 
<<<<<<< HEAD
  output$dtable <-  renderDataTable({
      dResult()
  })
  
  output$dBdownload <- downloadHandler(
      
      filename = function() { paste("drugbank_result.txt") },
      
      content = function(file) {
          write.table(dResult(), file, row.names=FALSE, quote=FALSE, sep="\t")
          
      })
 #addPopover(session, "Result", "Predicted Results", placement = "top",content = paste0("Shows the predicted results in a data table format"), trigger = 'click')
 #addPopover(session, "prop_table", "Network Properties", placement = "top",content = paste0("Represents the current selected network properties"), trigger = 'click')
 #addPopover(session, "advTable", "Network plot", placement = "top",content = paste0("This panel shows the predicted values"), trigger = 'click')
=======
 addPopover(session, "Result", "Predicted Results", placement = "top",content = paste0("Shows the predicted results in a data table format"), trigger = 'click')
 addPopover(session, "prop_table", "Network Properties", placement = "top",content = paste0("Represents the current selected network properties"), trigger = 'click')
 addPopover(session, "advTable", "Network plot", placement = "top",content = paste0("This panel shows the predicted values"), trigger = 'click')
 
>>>>>>> origin/master
})
