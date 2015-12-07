library(shinyBS)
library(shinythemes)
library(shiny)
library(shinysky)
library(shinyjs)
library(networkD3)
source("help.R")
library(visNetwork)
library(shinyGridster)
shinyUI(navbarPage(theme =shinytheme("spacelab"),img(src = "netpredicter.png", height = 32, width = 35) ,id ="bar",
                   tabPanel("home",
                            
                            div(style = " text-align: center;font-family: 'times'",
                                h2(" NetPredictor brings to you the power of finding missing links in Chem-Biological network .."),
                                br(),
                                br(),
                                img(src = "netpredicter.png", height = 150, width = 150,align = "success")
    
                                ),                   
                            br(),
                            br(),
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
                   tabPanel("Start Prediction",
                               
                            tags$head(
                                 tags$style(HTML("
                                  h1 {
                                    font-family: 'Verdana', cursive;
                                    line-height: 1.1;
                                    font-size : 20px;
                                    color: #4863A0
                                  }"))),
                            
                            headerPanel("Apply NetPredictor to your Data"),
                            
                            br(),br(),
#                             fluidRow(      
#                                 
#                                 column(4,wellPanel(
                            sidebarPanel(width=3,
                                # From http://stackoverflow.com/questions/19777515/r-shiny-mainpanel-display-style-and-font
                                HTML('<style type="text/css">
                                        .well { background-color: #ffffff; }
                                 </style>'),
                                       h4("load your data (or select Example)"),
                                       radioButtons(inputId="data_input_type", 
                                                    label="",
                                                    choices=c("Custom Data"="custom", "Example Data"="example"),
                                                    selected="", inline=FALSE),
                                       
                                      bsTooltip("data_input_type", "Choose Either example sets or input custom sets of data",
                                          "right", options = list(container = "body")),
                                
                                       ## Example data
                                       conditionalPanel(condition = "input.data_input_type == 'example'",
                                                    
                                                            selectInput("datasets", "Example datasets", choices = c("Enzyme","GPCR","Ion Channel","Nuclear Receptor"),selected = "",multiple = FALSE),
                                                        h4("Select Network Alogirthms"),
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
                                       actionButton('start', label='Start Prediction',
                                                     class="btn btn-primary"),
                                     render_helpfile("Data Input", "mds/import.md")
                                       ),mainPanel( 
                                           tabsetPanel(id='datatabs',
                                                       
                                                       tabPanel("Network Properties",
                                                                h3(textOutput("Data Summary", container = span)),
                                                                gridster(width = 250, height = 250,gridsterItem(row =1, col = 1, sizex = 1, sizey = 1, class = 'widget',
                                                                             uiOutput("prop_table")
                                                                ),gridsterItem(row =1, col = 2, sizex = 1, sizey = 1, class = 'widget',
                                                                               uiOutput("countDrugs")
                                                                ),gridsterItem(row =1, col = 3, sizex = 1, sizey = 1, class = 'widget',
                                                                               uiOutput("countProteins"))
                                                                )
                                                                ),
                                                       tabPanel("Prediction Results",
            
                                                                h4(textOutput("Prediction Results",container = span)),
            
                                                                dataTableOutput("Result"),
                                                                downloadButton("downloadResult", "Download results as csv file")),
                                                       tabPanel("networkplot",
                                                       h4("Network"),
                                                       visNetworkOutput("networkplot",height="700px"),
                                                       downloadButton("graphResult","Download Graph GML file"))
                                       )
                                       
                                       )),

                  ## Statistical Analysis TAB

                   navbarMenu(id='predtab',title="Advanced Analysis",
                   tabPanel("Statistical Analysis",value = "aa",
                        
                            
                            tags$head(
                                tags$style(HTML("
                                                h1 {
                                                font-family: 'Verdana', cursive;
                                                line-height: 1.1;
                                                font-size : 20px;
                                                color: #4863A0
                                                }"))),
                            
                            headerPanel("Perform Advanced Analysis on your Network"),
                            br(),br(),
                            sidebarPanel(width = 3,
                                # From http://stackoverflow.com/questions/19777515/r-shiny-mainpanel-display-style-and-font
                                HTML('<style type="text/css">
                                     .well { background-color: #ffffff; }
                                     </style>'),
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
                                             class="btn btn-primary"),
                                render_helpfile("Advanced Analysis", "mds/analysis.md")
                            ),mainPanel(tabPanel('Statistical Analysis',
                                                     h3(textOutput("Analysis")),
                                                     dataTableOutput("advTable")))
                            

                            
                            ),
                   
                   ## Permutation analysis Navbar menu tab
                   tabPanel("Permutation testing",
                            tags$head(
                                tags$style(HTML("
                                                h1 {
                                                font-family: 'Verdana', cursive;
                                                line-height: 1.1;
                                                font-size : 20px;
                                                color: #4863A0
                                                }"))),
                   headerPanel("Perform Random Permutation test on your network"),br(),
                   sidebarPanel(width = 3,
                                # From http://stackoverflow.com/questions/19777515/r-shiny-mainpanel-display-style-and-font
                                HTML('<style type="text/css">
                                     .well { background-color: #ffffff; }
                                     </style>'),
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
                                             class="btn btn-primary"),
                                render_helpfile("Significance Analysis", "mds/significance.md")
                   ),mainPanel(tabPanel('Permutation Analysis',
                                        h3(textOutput("Permutation Analysis")),
                                        dataTableOutput("sigTable"),
                                        downloadButton("downloadSig", "Download results as csv file")))
                   
                   
                   
                   )    
                   ),
             tabPanel("About",
                      HTML(markdown::markdownToHTML("mds/about.md", fragment.only=TRUE, options=c("")))
             )
                      
              


)
)




