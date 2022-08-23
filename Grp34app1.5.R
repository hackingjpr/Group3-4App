source("./AppSourceFunctions.R")

ui <- shiny::fluidPage(
  
  #Loading bar
  useAttendant(),
  #Dashboard looks etc
  dashboardPage(dashboardHeader(title= "", dropdownMenuOutput("messageMenu")), 
                dashboardSidebar(
                  #logo
                  # (img(src='Free_Sample_By_Wix (5).jpg', align = "center")),
                  (img(src='Free_Sample_By_Wix%20(12).jpg', align = "center")),
                  h2("Step 1"),
                  radioButtons(
                    "expmeth",
                    "Expression or Methylation Data?",
                    c(Expression = "exp", Methylation = "meth")
                  ),
                  conditionalPanel(
                    condition = "input.expmeth == 'meth'",
                    fileInput(
                      "methFile",
                      "Upload Methylation files:",
                      multiple = TRUE,
                      accept = ".idat")),
                    conditionalPanel(
                    condition = "input.expmeth == 'exp'",
                    fileInput(
                      "expFile",
                      "Upload Expression files:",
                      multiple = TRUE,
                      accept = c(".rds", ".csv", ".txt"))),
                  conditionalPanel(
                    condition = "input.expmeth == 'exp'",
                    sliderInput(
                      "outlier",
                      "Select outlier etc",
                      min = 0.5,
                      max = 10,
                      value = 0.5,
                      step = 0.5,
                      round = 0.5,
                      animate = TRUE)),
                  #uploading idat
                  # h2("Step 2"),
                  # fileInput(
                  #   "idatFile",
                  #   "Upload idat files:",
                  #   multiple = TRUE,
                  #   accept = ".idat"),
                  #Step 2
                  # h2("Step 2"),
                  # selectInput("metagenes", "Select Metagene Set", c("MRT (ATRT & ECRT)", "ATRT", "ECRT")),
                  #set it away
                  h2("Step 2"),
                  conditionalPanel(
                    condition = "input.expmeth == 'exp'",
                  actionBttn(
                    inputId = "bttn1",
                    label = "Generate Group 3/4 Score, Expression",
                    color = "primary",
                    style = "fill"
                  )),
                  conditionalPanel(
                    condition = "input.expmeth == 'meth'",
                    actionBttn(
                      inputId = "bttn2",
                      label = "Generate Group 3/4 Score, Methylation",
                      color = "success",
                      style = "fill"
                    )),
                  
                  #Reset button
                  actionBttn(
                    inputId = "bttn3",
                    label = "Reset",
                    color = "danger",
                    style = "fill",
                    size = "sm"
                  ),
                  #Loading Bar
                  attendantBar("progress-bar"),
                  hr(),
                  tags$footer( HTML(
                    paste(
                      "Author: James Hacking,", "<br/>",  
                      "Date Created: 26-07-2022,", "<br/>", 
                      "Copyright (c) James Hacking, 2022,", "<br/>",
                      h5(a("Email: james.hacking@ncl.ac.uk", href="mailto:james.hacking@ncl.ac.uk"))
                      
                    )))
                ), dashboardBody((
                  tabsetPanel( id = "tabs",
                               #Information Tab
                               tabPanel(
                                 "Info",
                                 (fluidRow(
                                   includeMarkdown("./introduction.md")
                                 )
                                 )
                               ),
                               #Tutorial Tab
                               tabPanel("Tutorial",
                                        includeMarkdown("./Tutorial/tutorial.md")
                               ),
                               #Results Tab
                               tabPanel("Results Table",
                                        (fluidRow(
                                          box(
                                            width = 12,
                                            title = "Group 3/4 score",
                                            status = "info",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            DTOutput('Mval')
                                          ))),
                                        # (fluidRow(
                                        #   box(
                                        #     title = "Risk plot", 
                                        #     status = "warning", 
                                        #     solidHeader = TRUE,
                                        #     collapsible = TRUE,
                                        #     plotOutput("figure")),
                                        #   box(
                                        #     title = "Selections",
                                        #     status = "success",
                                        #     solidHeader = TRUE,
                                        #     collapsible = TRUE,
                                        #     h3("Metagene Set:"),
                                        #     textOutput("metagenechoice"),
                                        #     h3("Sample Selected:"),
                                        #     textOutput("sample"),
                                        #     h3("Patient's Risk Percentile:"),
                                        #     textOutput("percentages")
                                        #   ))),
                                         (fluidRow(
                                           box(
                                             width = 12,
                                             title = "Summary",
                                             status = "primary",
                                             solidHeader = TRUE,
                                             collapsible = TRUE,
                                             includeMarkdown("./ResultsSummary.md")
                                           )
                                         ))
                                        ),
                               
                               #Download Tab
                               tabPanel("Download",
                                        textInput("filename", "Please insert desired filename", "Group 3/4 Score"),
                                        radioButtons(inputId = "download", label = "Select file type", choices = c("csv", "pdf")),
                                        downloadButton("down", "Download the results"),
                               )
                  )
                )),
                #Font Selection            
                # tags$head(tags$style(HTML('* {font-family: "Courier New"};')))),
  tags$head(tags$style(HTML('* {font-family: "Calibri"};')))),

  
)

server <- function(session, input, output) {
  # Allow larger files to be uploaded
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  #Hide tabs until graph made
  hideTab(inputId = "tabs", target = "Results Table")
  hideTab(inputId = "tabs", target = "Download")
  
  #Loading bar
  att <- Attendant$new("progress-bar")
  unlink("./temp/", recursive = T)
  if(!dir.exists("./temp")){
    dir.create("./temp")
  }
  #fix the overwriting issue
  fs::path("./temp/",createRandString()) -> tempDIR
  dir.create(tempDIR)
  
  Mvals <-
    observeEvent(input$bttn1, {
      req(input$expFile)
      att$set(0, text = "Loading") #Start at 10%
      att$auto(ms = 1600, value = 0.0048) # automatically increment
      message("Loading Bar")
      
      
      
      require(R.utils)
      
      input$expFile$datapath -> in.files
      message(in.files)
      paste0(tempDIR,"/", input$expFile$name) -> out.files
      message(out.files)
      message("TempDir")
      
      copied <- file.copy(in.files, out.files)
      message("copy")
      
      in.files <- readRDS(in.files)
      # in.files <- tpms.mat
      # message(head(in.files))
      
      # apply(g3g4,2,function(x){(1 / (1 + exp(-x)))}) -> logistic.g3g4
      # message("logistic.g3g4")
      # 
      # apply(logistic.g3g4,1,function(x){x[2]/(x[1]+x[2])}) -> logistic.g3g4.score
      # message("logistic.g3g4.score")
      # 
      # scaling.function(logistic.g3g4.score) -> logistic.g3g4.score
      # message("logistic.g3g4.score scaled")
      # 
      # # annotate nmb.mat to change into hgnc symbols
      # annotate.HTseq.IDs(rownames(nmb.mat)) -> annotation
      # message("annotation")
      # 
      # nmb.mat[-which(annotation$hgnc_symbol==""|is.na(annotation$hgnc_symbol)),] -> nmb.mat
      # message("nmb.mat")
      # 
      # annotation$hgnc_symbol[-which(annotation$hgnc_symbol==""|is.na(annotation$hgnc_symbol))] -> rownames(nmb.mat)
      # 
      # 
      # 
      # message("1")
      # 
      # 
      # ## pre-filter datasets
      # idx <- param.filter(2^nmb.mat)
      # nmb.mat <- 2^nmb.mat[idx, ]
      # message("2")
      # 
      # ## pre-projection normalisation
      # nmb.mat <- prep.data(nmb.mat)
      # message("3")
      
      nmb.mat <- nmb.mat.prepped
      # tpms.mat <- input$expFile$datapath
      
      
      # in.files <- as.matrix(in.files)
      print(object.size(in.files))  
      ## interset common genes / probes
      tpms.mat <- match.select(nmb.mat, in.files)
      # message(head(tpms.mat))
      message("4")
      

      ## NEED PACKAGE NMF
      init <- NMF::nmfModel(4, 
                            nmb.mat, 
                            W = 0.5, 
                            H = t(avg.h.val))
      message("5")
      
      ## generate NMF seeded with model
      nmf.res <- NMF::nmf(nmb.mat, 
                          4, 
                          seed = init, 
                          nrun = 8, 
                          .pbackend = 20)
      message("6")
      
      
      ## project using pseudo-inverse & post-projection normalise
      # project back onto the same dataset
      rnaseq.H <- project.NMF(input.array = nmb.mat, 
                              nmf.result = nmf.res)
      message("7")
      test <- as.matrix(tpms.mat)
      message(object.size(test))
      message(dim(test))
      # project onto fresh dataset
      tpms.H <- project.NMF(input.array = test, 
                            nmf.result = nmf.res)
      message("8")
      # tpms.H <- project.NMF(input.array = tpms.mat, 
      #                       nmf.result = nmf.res)
      
      ### define new g3g4 score for projection back onto the original data
      t(rnaseq.H[c(3,1),]) -> g3g4.rnaseq
      message("9")
      apply(g3g4.rnaseq,2,function(x){(1 / (1 + exp(-x)))}) -> logistic.g3g4.rnaseq
      message("10")
      apply(logistic.g3g4.rnaseq,1,function(x){x[2]/(x[1]+x[2])}) -> logistic.g3g4.rnaseq.score
      message("11")
      scaling.function(logistic.g3g4.rnaseq.score) -> logistic.g3g4.rnaseq.score
      message("12")
      
      ## join and plot the two together (original g3g4 score and g3g4 score projected back onto the same data) kind of a control that it is working
      ## NEED TIDYFT

      df <- inner_join(data.frame(logistic.g3g4.score,ids = names(logistic.g3g4.score)),
                       data.frame(logistic.g3g4.rnaseq.score, ids = names(logistic.g3g4.rnaseq.score)))
      message("13")
      
      t(tpms.H[c(3,1),]) -> g3g4.tpms
      message("14")
      apply(g3g4.tpms,2,function(x){(1 / (1 + exp(-x)))}) -> logistic.g3g4.tpms
      message("15")
      apply(logistic.g3g4.tpms,1,function(x){x[2]/(x[1]+x[2])}) -> logistic.g3g4.tpms.score
      message("16")
      
      # some times helpful to remove outliers prior to scaling
      outlier.idx <- c(head(order(logistic.g3g4.tpms.score), round((input$outlier)*(length(logistic.g3g4.tpms.score)/100))),
                       tail(order(logistic.g3g4.tpms.score), round((input$outlier)*(length(logistic.g3g4.tpms.score)/100)))
      )
      message("17")
      
      # scale to create final score
      scaling.function(logistic.g3g4.tpms.score[-outlier.idx]) -> logistic.g3g4.tpms.score
      round(logistic.g3g4.tpms.score, digits = 3) -> logistic.g3g4.tpms.score
      as.data.frame(logistic.g3g4.tpms.score) -> logistic.g3g4.tpms.score
      message("18")
      output$Mval <- renderDT (({logistic.g3g4.tpms.score
      }),
      options = list(
        pageLength = 10, 
        processing=FALSE),
      selection = "single"
      )
      
      showTab(inputId = "tabs", target = "Results Table", select = TRUE)
      showTab(inputId = "tabs", target = "Download")
    })
    
    
    observeEvent(input$bttn2, {
      att$set(10, text = "Loading") #Start at 10% 
      att$auto(ms = 1600, value = 0.01) # automatically increment
      
      #need idat file to run
      req(input$methFile)
      require(R.utils)
      
      
      input$methFile$datapath -> in.files
      message(in.files)
      paste0(tempDIR,"/", input$methFile$name) -> out.files
      message(out.files)
      
      copied <- file.copy(in.files, out.files)
      
      
      temp.base <- get_basenames(tempDIR)
      
      cat("Timing start\n")
      ptm <- proc.time()
      
      temp.processed <- process_idats(temp.base)
      on.exit({
        att$done()
      })
      

      output$Mval <- renderDT (
        ({
          # Obtain MValues
          beta2m(temp.processed$betas) -> M.values
          
          ### Round results to 3 figures
          # round(M.values, digits = 3) -> M.values
          
          metagene <- round(predict(g3.g4.cont.rfe, t(M.values)[,predictors(g3.g4.cont.rfe)]), digits = 3)
          metagene.df <- data.frame(g3g4.score = metagene)
          metagene.df
        }),
      options = list(
        pageLength = 10, 
        processing=FALSE),
      selection = "single"
      )

      att$done(text = "Complete")
      
      
      # output$time <- renderText({proc.time() - ptm})
      rowSelect <- reactive({input$Mval_rows_selected})
      message("row select done")
      #####
      # output$figure <-
      #   renderPlot({
      #     beta2m(temp.processed$betas) -> M.values
      #     pred.cont.rand.for <- predict(g3.g4.cont.rfe, t(M.values)[,predictors(g3.g4.cont.rfe)])
      #     round(pred.cont.rand.for, digits = 3) -> pred.cont.rand.for
      #     pred.cont.rand.for
      #       
      #       # print(figure.input)
      #     
      #     # generate_figure_highlight_mrt(
      #     #   figure.input
      #     #   ,1)
      #     
      #     # generate_figure_highlight_mrt(
      #     #   figure.input
      #     #   ,input$Mval_row_last_clicked)
      #     # 
      #     figure.output.mrt <- (
      #       (if (input$metagenes == "MRT (ATRT & ECRT)") {
      #       # figureFile <- "./mrt54.dist.rds"
      #       generate_figure_highlight_mrt(
      #         figure.input
      #         ,input$Mval_row_last_clicked)
      #     }))
      #     figure.output.atrt <- (
      #       (if(input$metagenes == "ATRT") {
      #       # figureFile <- "./atrt8.dist.rds"
      #       generate_figure_highlight_atrt(
      #         figure.input
      #         ,input$Mval_row_last_clicked)
      #     }))
      #     figure.output.ecrt <-( 
      #       (if(input$metagenes == "ECRT") {
      #       # figureFile <- "./ecrt20.dist.rds"
      #       generate_figure_highlight_ecrt(
      #         figure.input
      #         ,input$Mval_row_last_clicked)
      #     }))
      #     # }
      #     
      #     if (input$metagenes == "MRT (ATRT & ECRT)") {
      #       figure.output.mrt -> figure.output
      #     }
      #     if (input$metagenes == "ATRT") {
      #       figure.output.atrt -> figure.output
      #     }
      #     if (input$metagenes == "ECRT") {
      #       figure.output.ecrt -> figure.output
      #     }
      #     
      #     
      #     figure.output
      #     # https://rstudio.github.io/DT/shiny.html
      #     
      #     
        # })
      # output$percentages <- renderText (
      #   {
      #   test.res.mrt <- extract.metagene(
      #     as.character(ALL[[1]]$genes),
      #     as.numeric(ALL[[1]]$weights),
      #     beta2m(temp.processed$betas),
      #     as.numeric(ALL[[2]])
      #   )
      #   
      #   round(test.res.mrt, digits = 3) -> test.res.mrt 
      #   
      #   
      #   test.res.atrt <- extract.metagene(
      #     as.character(ATRT[[1]]$genes),
      #     as.numeric(ATRT[[1]]$weights),
      #     beta2m(temp.processed$betas),
      #     as.numeric(ATRT[[2]])
      #   )
      #   
      #   round(test.res.atrt, digits = 3) -> test.res.atrt
      #   
      #   
      #   test.res.ecrt <- extract.metagene(
      #     as.character(ECRT[[1]]$genes),
      #     as.numeric(ECRT[[1]]$weights),
      #     beta2m(temp.processed$betas),
      #     as.numeric(ECRT[[2]])
      #   )
      #   
      #   round(test.res.ecrt, digits = 3) -> test.res.ecrt
      #   
      #   if (input$metagenes == "MRT (ATRT & ECRT)") {
      #     test.res.mrt$Risk_Value -> figure.input
      #     names(figure.input) <- rownames(test.res.mrt)
      #   }
      #   if (input$metagenes == "ATRT") {
      #     test.res.atrt$Risk_Value -> figure.input
      #     names(figure.input) <- rownames(test.res.atrt)
      #   }
      #   if (input$metagenes == "ECRT") {
      #     test.res.ecrt$Risk_Value -> figure.input
      #     names(figure.input) <- rownames(test.res.ecrt)
      #   }
      #   
      #   # print(figure.input)
      #   
      #   
      #   figure.output.mrt <- (
      #     (if (input$metagenes == "MRT (ATRT & ECRT)") {
      #       # figureFile <- "./mrt54.dist.rds"
      #       generate_figure_percentage_mrt(
      #         figure.input
      #         ,input$Mval_row_last_clicked)
      #     }))
      #   figure.output.atrt <- (
      #     (if(input$metagenes == "ATRT") {
      #       # figureFile <- "./atrt8.dist.rds"
      #       generate_figure_percentage_atrt(
      #         figure.input
      #         ,input$Mval_row_last_clicked)
      #     }))
      #   figure.output.ecrt <-( 
      #     (if(input$metagenes == "ECRT") {
      #       # figureFile <- "./ecrt20.dist.rds"
      #       generate_figure_percentage_ecrt(
      #         figure.input
      #         ,input$Mval_row_last_clicked)
      #     }))
      #   # }
      #   
      #   if (input$metagenes == "MRT (ATRT & ECRT)") {
      #     figure.output.mrt -> figure.output.percentage
      #   }
      #   if (input$metagenes == "ATRT") {
      #     figure.output.atrt -> figure.output.percentage
      #   }
      #   if (input$metagenes == "ECRT") {
      #     figure.output.ecrt -> figure.output.percentage
      #   }
      #   
      #   
      #   figure.output.percentage
      #   
      #   # input <- "./mrt54.dist.rds",
      #   # generate_figure_percentage_mrt(
      #   #   figure.input,
      #   #   input$Mval_row_last_clicked)
      #   # ,input$Mval_row_last_clicked)
      #   # "hello",
      #   
      #    # if (input$metagenes == "MRT (ATRT & ECRT)") {
      #    #   # figureFile <- "./mrt54.dist.rds"
      #    #   generate_figure_percentage_mrt(
      #    #     figure.input
      #    #     ,input$Mval_row_last_clicked)
      #    # }
      #    # if (input$metagenes == "ATRT") {
      #    #   # figureFile <- "./atrt8.dist.rds"
      #    #   generate_figure_percentage_atrt(
      #    #     figure.input
      #    #     ,input$Mval_row_last_clicked)
      #    # }
      #    # if (input$metagenes == "ECRT") {
      #    #   # figureFile <- "./ecrt20.dist.rds"
      #    #   generate_figure_percentage_ecrt(
      #    #     figure.input
      #    #     ,input$Mval_row_last_clicked)
      #    # }
      # })
#####
      # output$metagenechoice <- renderText({input$metagenes})
      
      showTab(inputId = "tabs", target = "Results Table", select = TRUE)
      showTab(inputId = "tabs", target = "Download")
      
      #For displaying currently selected sample
      
      output$sample <- renderText(input$Mval_cell_clicked$value) 
      
      output$down <- downloadHandler(
        filename = function() {
          paste(input$metagenes,Sys.time(), input$download, sep=".")
        },
        content ={ 
          function(file) {
            if (input$download == "csv")
            write.csv(
              (
                ({
                  beta2m(temp.processed$betas) -> M.values
                  
                  ### Round results to 3 figures
                  # round(M.values, digits = 3) -> M.values
                  
                  metagene <- round(predict(g3.g4.cont.rfe, t(M.values)[,predictors(g3.g4.cont.rfe)]), digits = 3)
                  metagene.df <- data.frame(g3g4.score = metagene)
                  metagene.df
                  }) 
                
              ), file, row.names = TRUE)
            ####### if not csv then write a pdf
            else
              pdf(file,
                  width = 14,
                  title = paste(input$metagenes, "output")
                  )
            
            grid.table(
              
              ({
                beta2m(temp.processed$betas) -> M.values
                
                ### Round results to 3 figures
                # round(M.values, digits = 3) -> M.values
                
                metagene <- round(predict(g3.g4.cont.rfe, t(M.values)[,predictors(g3.g4.cont.rfe)]), digits = 3)
                metagene.df <- data.frame(g3g4.score = metagene)
                metagene.df
                })
            )
            #####
            # plot({
            #   
            #   
            #     
            #     test.res.mrt <- extract.metagene(
            #       as.character(ALL[[1]]$genes),
            #       as.numeric(ALL[[1]]$weights),
            #       beta2m(temp.processed$betas),
            #       as.numeric(ALL[[2]])
            #     )
            #     
            #     round(test.res.mrt, digits = 3) -> test.res.mrt 
            #     
            #     
            #     test.res.atrt <- extract.metagene(
            #       as.character(ATRT[[1]]$genes),
            #       as.numeric(ATRT[[1]]$weights),
            #       beta2m(temp.processed$betas),
            #       as.numeric(ATRT[[2]])
            #     )
            #     
            #     round(test.res.atrt, digits = 3) -> test.res.atrt
            #     
            #     
            #     test.res.ecrt <- extract.metagene(
            #       as.character(ECRT[[1]]$genes),
            #       as.numeric(ECRT[[1]]$weights),
            #       beta2m(temp.processed$betas),
            #       as.numeric(ECRT[[2]])
            #     )
            #     
            #     round(test.res.ecrt, digits = 3) -> test.res.ecrt
            #     
            #     if (input$metagenes == "MRT (ATRT & ECRT)") {
            #       test.res.mrt$Risk_Value -> figure.input
            #       names(figure.input) <- rownames(test.res.mrt)
            #     }
            #     if (input$metagenes == "ATRT") {
            #       test.res.atrt$Risk_Value -> figure.input
            #       names(figure.input) <- rownames(test.res.atrt)
            #     }
            #     if (input$metagenes == "ECRT") {
            #       test.res.ecrt$Risk_Value -> figure.input
            #       names(figure.input) <- rownames(test.res.ecrt)
            #     }
            #     
            #     # print(figure.input)
            #     
            #     # generate_figure_highlight_mrt(
            #     #   figure.input
            #     #   ,1)
            #     
            #     # generate_figure_highlight_mrt(
            #     #   figure.input
            #     #   ,input$Mval_row_last_clicked)
            #     # 
            #     figure.output.mrt <- (
            #       (if (input$metagenes == "MRT (ATRT & ECRT)") {
            #         # figureFile <- "./mrt54.dist.rds"
            #         generate_figure_highlight_mrt(
            #           figure.input
            #           ,NA)
            #       }))
            #     figure.output.atrt <- (
            #       (if(input$metagenes == "ATRT") {
            #         # figureFile <- "./atrt8.dist.rds"
            #         generate_figure_highlight_atrt(
            #           figure.input
            #           ,NA)
            #       }))
            #     figure.output.ecrt <-( 
            #       (if(input$metagenes == "ECRT") {
            #         # figureFile <- "./ecrt20.dist.rds"
            #         generate_figure_highlight_ecrt(
            #           figure.input
            #           ,NA)
            #       }))
            #     # }
            #     
            #     if (input$metagenes == "MRT (ATRT & ECRT)") {
            #       figure.output.mrt -> figure.output
            #     }
            #     if (input$metagenes == "ATRT") {
            #       figure.output.atrt -> figure.output
            #     }
            #     if (input$metagenes == "ECRT") {
            #       figure.output.ecrt -> figure.output
            #     }
            #     
            #     
            #     figure.output
            #     
            #     
            #   })####
              dev.off()
            
          }
          # if (input$download == "csv")####
          #   write.csv(output$Mval)
          # else
          #   pdf(output$Mval,
          #       width = 14
          #   )
          # dev.off()
        })
      
      
    })
  #Reset session and delete tempDIR
  observeEvent(input$bttn3, {
    unlink(tempDIR, recursive = T)
    session$reload()
    return()
    print("session reload not working")
  })
  renderText(output$metagenes <- input$metagenes)
  #####
  # output$down <- downloadHandler(
  #   filename = function() {
  #     paste(input$filename,Sys.time(), ".csv", sep="_")
  #   },
  #   content ={ 
  #     function(file) {
  #       write.csv(
  #         (
  #           ({
  #             test.res.mrt <- extract.metagene(
  #               as.character(ALL[[1]]$genes),
  #               as.numeric(ALL[[1]]$weights),
  #               beta2m(temp.processed$betas),
  #               as.numeric(ALL[[2]])
  #             )
  #             
  #             round(test.res.mrt, digits = 3) -> test.res.mrt
  #             
  #             
  #             test.res.atrt <- extract.metagene(
  #               as.character(ATRT[[1]]$genes),
  #               as.numeric(ATRT[[1]]$weights),
  #               beta2m(temp.processed$betas),
  #               as.numeric(ATRT[[2]])
  #             )
  #             
  #             round(test.res.atrt, digits = 3) -> test.res.atrt
  #             
  #             
  #             test.res.ecrt <- extract.metagene(
  #               as.character(ECRT[[1]]$genes),
  #               as.numeric(ECRT[[1]]$weights),
  #               beta2m(temp.processed$betas),
  #               as.numeric(ECRT[[2]])
  #             )
  #             
  #             round(test.res.ecrt, digits = 3) -> test.res.ecrt
  #             
  #             if (input$metagenes == "MRT (ATRT & ECRT)") {
  #               test.res.mrt -> test.res
  #             }
  #             if (input$metagenes == "ATRT") {
  #               test.res.atrt -> test.res
  #             }
  #             if (input$metagenes == "ECRT") {
  #               test.res.ecrt -> test.res
  #             }
  #             
  #             test.res}) 
  #           
  #         ), file, row.names = TRUE)
  #     }
  #     # if (input$download == "csv")
  #     #   write.csv(output$Mval)
  #     # else
  #     #   pdf(output$Mval,
  #     #       width = 14
  #     #   )
  #     # dev.off()
  #   })
  # 
  
  
  
  
  
}

shinyApp(server = server, ui = ui)