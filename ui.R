# Copyright (C) 2016 Abhik Seal <abhik1368@gmail.com>
# This program CC BY-NC-SA 3.0 license.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#################################################################################
## Load the libraries
library(shinyBS)
library(shinythemes)
library(shiny)
library(shinysky)
library(shinyjs)
library(rCharts)
source("help.R")
library(visNetwork)
library(shinyLP)
library(V8)
library(biomaRt)
library(bsplus)
library(DT)
options(shiny.trace = TRUE)
options(shiny.reactlog = TRUE)
options(shiny.fullstacktrace = TRUE)
library(shinyWidgets)
#library(shinyGridster)
#################################################################################

inputTextarea <- function(inputId, value="", nrows, ncols) {
    tagList(
        singleton(tags$head(tags$script(src = "textarea.js"))),
        tags$textarea(id = inputId,
                      class = "inputtextarea",
                      rows = nrows,
                      cols = ncols,
                      as.character(value))
    )
}
jsCode <- "shinyjs.browseURL = function(URL){window.open(URL,''); ;}"
jscode <- "shinyjs.closeWindow = function() { window.close(); }"

navbarPageWithText <- function(..., text) {
  navbar <- navbarPage(...)
  textEl <- tags$p(class = "navbar-text pull-right", text)
  navbar[[3]][[1]]$children[[1]] <- htmltools::tagAppendChild(
    navbar[[3]][[1]]$children[[1]], textEl)
  navbar
}



  
shinyUI(navbarPage(theme = shinytheme("flatly"),img(src = "netpredicter.png", height = 32, width = 35) ,id ="bar",
                   tabPanel(icon("home",lib = "glyphicon"),
                            #includeCSS("www/materialize.min.css"),
                            div(style = " text-align: center;font-family: 'Roboto'",
                                h2(" NetPredictor brings to you the power of finding missing links in Chem-Biological network .."),
                                br(),
                                br(),
                                img(src = "netpredicter.png", height = 150, width = 150,align = "success")

                                ),
                            br(),
                            # br(),
                            br(),
                            div(style = "text-align: center;",
                                bsButton("gofind", label = "Start Prediction",  style = "default",
                                         size = "large", disabled = FALSE
                                )),
                            br(),br(),br(),br(),br(),br(),
                            div(style= "text-align: center;",
                                fluidRow(
                                    column(4,
                                           img(src = "upload.png", height = 90, width = 100,align = "center"),
                                           h3("1. Load your data into NetPredicter ")),
                                    column(4,
                                           img(src = "network.png", height = 90, width = 85,align = "center"),
                                           h3("2. Predict for missing links in a Network")),

                                    column(4,
                                           img(src = "discover.png", height = 90, width = 100,align = "center"),

                                           h3("3. Discover new relations in a network")
                                    )),
                                br(),
                                br()
                            )
                   ),
                   tabPanel("Start Prediction",icon=icon("dot-circle-o"),

                            # tags$head(
                            #      tags$style(HTML("
                            #       h1 {
                            #         font-family: 'Verdana', cursive;
                            #         line-height: 1.1;
                            #         font-size : 20px;
                            #         color: #4863A0
                            #       }"))),

                            headerPanel("Apply NetPredictor to your Data"),

                            br(),br(),
#                             fluidRow(
#
#                                 column(4,wellPanel(
                            sidebarPanel(width=3,
                                         HTML('<style type="text/css">
                                              .well { background-color: #ffffff; }
                                              </style>'),
                                conditionalPanel(condition="input.conditionedPanels==1",
                                       h4(icon("upload",lib = "glyphicon"),"load your data (or select Example)"),
                                       radioButtons(inputId="data_input_type",
                                                    label="",
                                                    choices=c("Custom Data"="custom", "Example Data"="example"),
                                                    selected="", inline=FALSE),

                                      bsTooltip("data_input_type", "Choose Either example sets or input custom sets of data",
                                          "right", options = list(container = "body")),

                                       ## Example data
                                       conditionalPanel(condition = "input.data_input_type == 'example'",

                                                            selectInput("datasets", "Example datasets", choices = c("Enzyme","GPCR","Ion Channel","Nuclear Receptor"),selected = "",multiple = FALSE),
                                                        h4(icon("magic",lib = "font-awesome"),"Select Network Alogirthms"),
                                                        radioButtons(inputId="algorithm_typeI",
                                                                     label="",
                                                                     choices=c("HeatS"="heat","Network Based Inference"="nbi","Random Walk with Restart"="rwr","NetCombo"="netcombo"),
                                                                     selected="", inline=FALSE),
                                                        conditionalPanel(condition = "input.algorithm_typeI=='heat'",
                                                                         sliderInput("heatAlpha", "Alpha", value=0.5, min=0, max=1,width ='200px'),
                                                                         sliderInput("heatLambda", "Lambda", value=0.5, min=0, max=1,width = '200px')),
                                                        conditionalPanel(condition = "input.algorithm_typeI=='nbi'",
                                                                         sliderInput("nbiAlpha", "Alpha", value=0.5, min=0, max=1,width ='200px'),
                                                                         sliderInput("nbiLambda", "Lambda", value=0.5, min=0, max=1,width = '200px')),
                                                        conditionalPanel(condition = "input.algorithm_typeI=='rwr'",
                                                                         sliderInput("rwr_restart", "Restart", value=0.8, min=0, max=1,width ='200px')),
                                                        conditionalPanel(condition = "input.algorithm_typeI=='netcombo'",
                                                                         sliderInput("nc_restart", "Restart", value=0.8, min=0, max=1,width ='200px'),
                                                                         sliderInput("ncAlpha", "Alpha", value=0.5, min=0, max=1,width ='200px'),
                                                                         sliderInput("ncLambda", "Lambda", value=0.5, min=0, max=1,width ='200px'))),
                                       # Local files
                                       conditionalPanel(condition = "input.data_input_type == 'custom'",

                                                        fileInput('dt_file', label="Drug-Target Biparite Network",
                                                                  multiple=FALSE,
                                                                  accept=c('text/csv',
                                                                           'text/comma-separated-values',
                                                                           'text/tab-separated-values',
                                                                           '.csv',
                                                                           'text/plain',
                                                                           '.tsv')),
                                                        bsTooltip("dt_file", "Input CSV file with Binary Drug Target Interactions Matrix",
                                                                  "right", options = list(container = "body")),
                                                        fileInput('drug_file', label="Drug Similarity Matrix",
                                                                  multiple=FALSE,
                                                                  accept=c('text/csv',
                                                                           'text/comma-separated-values',
                                                                           'text/tab-separated-values',
                                                                           '.csv',
                                                                           'text/plain',
                                                                           '.tsv')),
                                                        bsTooltip("drug_file", "Input CSV file with Drug-Drug similarity Matrix",
                                                                  "right", options = list(container = "body")),
                                                        fileInput('target_file', label="Target Similarity Matrix",
                                                                  multiple=FALSE,
                                                                  accept=c('text/csv',
                                                                           'text/comma-separated-values',
                                                                           'text/tab-separated-values',
                                                                           '.csv',
                                                                           'text/plain',
                                                                           '.tsv')),
                                                        bsTooltip("target_file", "Input CSV file with Target-Target similarityMatrix",
                                                                  "right", options = list(container = "body")),
                                                        h4("Select Network Alogirthms"),
                                                        radioButtons(inputId="algorithm_typeII",
                                                                     label="",
                                                                     choices=c("HeatS"="heat","Network Based Inference"="nbi","Random Walk with Restart"="rwr","NetCombo"="netcombo"),
                                                                     selected="nbi", inline=FALSE),
                                                        conditionalPanel(condition = "input.algorithm_typeII=='heat'",
                                                                         sliderInput("cheatAlpha", "Alpha", value=0.5, min=0, max=1,width ='200px'),
                                                                         sliderInput("cheatLambda", "Lambda", value=0.5, min=0, max=1,width = '200px')),
                                                        conditionalPanel(condition = "input.algorithm_typeII=='nbi'",
                                                                         sliderInput("cnbiAlpha", "Alpha", value=0.5, min=0, max=1,width ='200px'),
                                                                         sliderInput("cnbiLambda", "Lambda", value=0.5, min=0, max=1,width = '200px')),
                                                        conditionalPanel(condition = "input.algorithm_typeII=='rwr'",
                                                                         sliderInput("crwr_restart", "Restart", value=0.8, min=0, max=1,width ='200px')),
                                                        conditionalPanel(condition = "input.algorithm_typeII=='netcombo'",
                                                                         sliderInput("cnc_restart", "Restart", value=0.8, min=0, max=1,width ='200px'),
                                                                         sliderInput("cncAlpha", "Alpha", value=0.5, min=0, max=1,width ='200px'),
                                                                         sliderInput("cncLambda", "Lambda", value=0.5, min=0, max=1,width ='200px'))


                                       ),
                                      busyIndicator("Calculation In progress",wait = 0),
                                       actionButton('start', label='Run Prediction',
                                                     class="btn btn-primary")
                                      ),
                                     conditionalPanel(condition="input.conditionedPanels==2",
                                                      h4("Get Predictive Metrics"),
                                                      numericInput('relinks',width = '200px',
                                                                   label = 'Choose Random links to be removed',
                                                                   min = 1, value = 25),
                                                      numericInput('freqT',width = '200px',
                                                                   label = 'Frequency of associations between Biparite Nodes',
                                                                   min = 1, value = 2),
                                                      radioButtons(inputId="predMetrics",
                                                                   label="",
                                                                   choices=c("Network Based Inference"="nbi","Random walk with restart"="rwr",
                                                                             "NetCombo"='nc'),
                                                                   selected="", inline=FALSE),
                                                      conditionalPanel(condition = "input.predMetrics=='nbi'",
                                                                       sliderInput("pdnbiAlpha", "Alpha", value=0.5, min=0, max=1,width ='200px'),
                                                                       sliderInput("pdnbiLambda", "Lambda", value=0.5, min=0, max=1,width = '200px')),
                                                      conditionalPanel(condition = "input.predMetrics=='rwr'",
                                                                       sliderInput("pdrwrRestart", "Restart", value=0.8, min=0, max=1,width ='200px')),
                                                      conditionalPanel(condition = "input.predMetrics=='nc'",
                                                                       sliderInput("pdcnc_restart", "Restart", value=0.8, min=0, max=1,width ='200px'),
                                                                       sliderInput("pdcncAlpha", "Alpha", value=0.5, min=0, max=1,width ='200px'),
                                                                       sliderInput("pdcncLambda", "Lambda", value=0.5, min=0, max=1,width ='200px')),
                                                      busyIndicator("Calculation In progress",wait = 0),
                                                      actionButton('submit', label='Submit',
                                                                   class="btn btn-primary")
                                                      # render_helpfile("Advanced Analysis", "mds/analysis.md")
                                     ),
                                   conditionalPanel(condition="input.conditionedPanels==3",
                                                    h4("Perform Permutations Analysis on your Network"),
                                                    numericInput('permute',width = '200px',
                                                                 label = 'Choose number of random permutations',
                                                                 min = 5, value = 10),
                                                    numericInput('sig',width = '150px',
                                                                 label = 'Keep Significant links of pvalue',
                                                                 min = 0.00000000001, value =0.05),
                                                    radioButtons(inputId="sigMetrics",
                                                                 label="",
                                                                 choices=c("Network Based Inference"="signbi","Random walk with restart"="sigrwr"),
                                                                 selected="", inline=FALSE),
                                                    conditionalPanel(condition = "input.sigMetrics=='signbi'",
                                                                     sliderInput("sgnbiAlpha", "Alpha", value=0.5, min=0, max=1,width ='200px'),
                                                                     sliderInput("sgnbiLambda", "Lambda", value=0.5, min=0, max=1,width = '200px')),
                                                    conditionalPanel(condition = "input.sigMetrics=='sigrwr'",
                                                                     sliderInput("sgrwrRestart", "Restart", value=0.8, min=0, max=1,width ='200px')),
                                                    busyIndicator("Calculation In progress",wait = 0),
                                                    actionButton('sigSubmit', label='Submit',
                                                                 class="btn btn-primary")

                                                    )

                                      ),mainPanel(
                                           tabsetPanel(id = "conditionedPanels",

                                                       tabPanel("Network Properties",br(),value = 1,
                                                                actionButton('netproperty', label='Calculate Properties',class="btn btn-primary"),
                                                                h3(textOutput("Data Summary", container = span)),
                                                                uiOutput("prop_table"),
                                                                tabsetPanel(
                                                                    id = 'dataset',tabPanel("Proteins Properties",DT::dataTableOutput("btwProteinsdt")),
                                                                    tabPanel("Drug Properties",DT::dataTableOutput("btwDrugsdt")))
                                                                ),
                                                                    #column(width = 5, chartOutput("countProteins","polycharts")),
                                                                    #column(width = 6, offset = 1, chartOutput("countDrugs","polycharts"))
                                                                #),
                                                      #  fluidRow(
                                                      #
                                                      #     #column(width = 5, chartOutput("btwProteins","polycharts")),
                                                      #     column(width = 6, offset = 1, chartOutput("btwDrugs","polycharts"))
                                                      # )),
#                                                                 showOutput("btwProteins","polycharts"),
#                                                                 showOutput("btwDrugs","polycharts")),

                                                       tabPanel("Network Modules",value = 1,

                                                                br(),
                                                                actionButton('mods', label='Calculate Modules',class="btn btn-primary"),
                                                                #render_helpfile("Network Modules", "mds/module.md"),
                                                                uiOutput('modules'),
                                                                dataTableOutput("data_table"),
                                                                actionButton('shownet', label='Show Network',class="btn btn-primary"),
                                                                div(style="width:1200px; height:225px;",bsCollapse(id = "module1", open = "modules",
                                                                                                                   bsCollapsePanel("modules",visNetworkOutput("moduleplot",height="550px"))))
                                                                ),
                                                       tabPanel("Prediction Results",value = 1,

                                                                h4(textOutput("Prediction Results",container = span)),

                                                                dataTableOutput("Result"),
                                                                downloadButton("downloadResult", "Download results as csv file")),
                                                       tabPanel("networkplot",value = 1,
                                                       h4("Network"),
                                                       visNetworkOutput("networkplot",height="725px"),
                                                       downloadButton("graphResult","Download Graph GML file")),
                                                       tabPanel('Statistical Analysis',value =2,
                                                                 h3(textOutput("Analysis")),
                                                                 dataTableOutput("advTable"),
                                                                 downloadButton("downloadadvr", "Download results as csv file")),
                                                        tabPanel('Permutation Analysis', value = 3,
                                                                 h3(textOutput("Permutation Analysis")),
                                                                 dataTableOutput("sigTable"),
                                                                 downloadButton("downloadSig", "Download results as csv file"))


                                       )

                                       )),

                  ## Statistical Analysis TAB
                tabPanel("PPI Network",
                         div(id = "demo", class = "collapse in", 
                             
                         sidebarPanel(width = 3,HTML('<style type="text/css">
                                              .well { background-color: #ffffff; }
                                                     </style>'),
                                      div(style="display:inline-block",h4(icon("search",lib = "glyphicon"),"Find Interactions")),
                                      radioButtons(inputId="data_search",
                                                   label="",
                                                   choices=c("Protein-Protein Interactions(PPI)"="ppi", "PPI Shortest Paths"="path","Subgraph Extract"="subgraph" ),
                                                   selected="", inline=FALSE),

                                      conditionalPanel(
                                        condition = "input.data_search == 'ppi'",
                                        selectizeInput("ppi", label= "HGNC Symbols", choices = NULL, selected = "", multiple = T,width='85%'),
                                        
                                        #uiOutput("location"),
                                        sliderInput("confidence", "Confidence:",
                                                    min = 0, max = 1, value = 0.99, step= 0.01),
                                        div(style="display:inline-block", radioButtons("dataid", "Select only:",  choices=c("show ppi", "show drugs"), selected ="show ppi"))),

                                        #uiOutput("ppipanel")),

                                       conditionalPanel(
                                        condition = "input.data_search == 'path'",
                                        radioButtons(inputId="weight",
                                                     label="",
                                                     choices=c("Weighted PPI", "Unweighted PPI"),
                                                     selected="Weighted PPI", inline=FALSE),

                                         textInput("from", label = "From :",
                                                              value = "",
                                                             width = "100px"),

                                        textInput("to", label = "To:",
                                                             value = "",
                                                              width = "100px"),
                                       numericInput("kpath", label = "Top-Paths :", min = 1,
                                                   max = 20, value = 1 , step = 1,
                                                   width = "150px"),
                                      checkboxInput('pathway', label = "Enrich Pathways (Reactome)", value = FALSE, width = NULL),
                                      uiOutput("conditionalCheck")),
                                      
                                      conditionalPanel(
                                        condition = "input.data_search == 'subgraph'",
                                        inputTextarea('proteinlist', '',8,20),br(),
                                        div(style="display:inline-block", radioGroupButtons(inputId = "ppi_id3", 
                                                                                            label = "Interactomes", choices = c("ConsensusPath", 
                                                                                                                                "String"), size="sm",status = "primary", selected = "ConsensusPath",
                                                                                            checkIcon = list(yes = icon("ok",lib = "glyphicon"))))),
                                      
                                     div(style="display:inline-block" ,actionButton('ppisearch', label='Search',class="btn btn-primary"))
                                     )
                         ),
                         mainPanel(conditionalPanel(condition = "input.data_search == 'ppi'",
                                                                   br(),
                                                                   bsAlert("alertppi"),
                                                    div(style="width:1200px; height:225px;",bsCollapse(id="collapse1",open="PPI Network",
                                                                                                       bsCollapsePanel("PPI Network",visNetworkOutput("ppinet",height="675px"))))),

                                  conditionalPanel(condition = "input.data_search == 'path'",
                                                                       bsAlert("alert"),
                                                   div(style="width:1200px; height:225px;",bsCollapse(id="collapse1",open="Shotest paths",bsCollapsePanel("Shotest paths",visNetworkOutput("pathnet",height="675px"))),
                                                   
                                                   bsCollapsePanel("pathways",DT::dataTableOutput("pathTab",width = "40%"), style = "font-size: 40%"))),
                                  conditionalPanel(condition = "input.data_search == 'subgraph'",
                                                   div(style="width:1200px; height:225px;",bsCollapse(id = "collapseExample", open = "PPI Subgraph",
                                                                                                      bsCollapsePanel("PPI Subgraph",visNetworkOutput("subnet",height = "700px")))))
                                             )

                         ),
                   tabPanel("Search Drugbank", icon= icon("search"),
                            tags$head(
                                tags$style(HTML("
                                                h1 {
                                                font-family: 'Verdana', cursive;
                                                line-height: 1.1;
                                                font-size : 16px;
                                                color: #4863A0
                                                }"))),
                            headerPanel("Select Drug names/drugbank IDs to search"),br(),
                            sidebarPanel(width = 3,
                                         # HTML('<style type="text/css">
                                         #      .well { background-color: #ffffff; }
                                         #      </style>'),
                                         radioButtons(inputId="search_type",
                                                      label = "",
                                                      choices=c("Search Drugs"="drugs", "Search Proteins"="proteins"),
                                                      selected="", inline=FALSE),
                                         conditionalPanel(condition = "input.search_type == 'drugs'",
                                         textInput("did", "Drugbank ID:", width = NULL)),
                                         conditionalPanel(condition = "input.search_type == 'proteins'",
                                                          textInput("pid", "Hugo Gene:", width = NULL)),
                                         busyIndicator("Search In progress",wait = 0),
#                                          radioButtons(inputId="algo_dtype",
#                                                       label="Algorithm Type",
#                                                       choices=c("NBI"="nbi", "RWR"="rwr"),
#                                                       selected="nbi", inline=FALSE),
                                         actionButton('dSearch', label='Submit',
                                                      class = "btn btn-primary")

                            ),mainPanel(tabPanel('Predicted Links',

                                   dataTableOutput("dtable"),
                                   downloadButton("dBdownload", "Download results as csv file")
                                   ))),
              tabPanel("Enrichment Analysis" , icon = icon("tasks"),
                       tags$head(
                           tags$style(HTML("
                                           h1 {
                                           font-family: 'Verdana', cursive;
                                           line-height: 1.1;
                                           font-size : 16px;
                                           color: #4863A0
                                           }"))),
                            headerPanel("Paste gene on a new line for enrichment analysis"),br(),
                            sidebarPanel(width = 3,
                                    # HTML('<style type="text/css">
                                    #      .well { background-color: #ffffff; }
                                    #      </style>'),
                                    inputTextarea('selectGene', '',8,15 ),
                                    radioButtons(inputId="gopath",
                                                 label="Select Ontology or Pathway or Disease",
                                                 choices=c("Gene Ontology"="go","Pathway Enrichment"="pathway","Disease Enrichment"="disease"),
                                                 selected="", inline=FALSE),
                                    conditionalPanel(condition = "input.gopath=='go'",
                                                     textInput("level", label = "Enter GO Level", value = 3)),
                                    busyIndicator("Search In progress",wait = 0),
                                    actionButton('genelist', label='Submit',
                                                 class="btn btn-primary")

                            ),mainPanel(tabPanel('Enrichment Results',

                                                 dataTableOutput("genePathway"),
                                                 downloadButton("Godownload", "Download results as csv file")
                            ))),
              tabPanel("About", icon= icon("info-circle"),
                      HTML(markdown::markdownToHTML("mds/about.md", fragment.only=TRUE, options=c("")))
             )


))






