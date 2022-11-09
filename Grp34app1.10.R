source("./AppSourceFunctions1.2.R")

library(shinythemes)

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
                      "OUTLIER SELECTION WHAT SHOULD I CALL THIS??",
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
                    style = "gradient"
                  )),
                  conditionalPanel(
                    condition = "input.expmeth == 'meth'",
                    actionBttn(
                      inputId = "bttn2",
                      label = "Generate Group 3/4 Score, Methylation",
                      color = "success",
                      style = "gradient"
                    )),
                  
                  #Reset button
                  actionBttn(
                    inputId = "bttn3",
                    label = "Reset",
                    color = "danger",
                    style = "minimal",
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
                      h5(a("Email: james.hacking@ncl.ac.uk", href="mailto:james.hacking@ncl.ac.uk")),
                      h4("Disclaimer : This app is designed exclusively for research purposes and is strictly not for diagnostic use.")
                      # shinythemes::themeSelector,()
                      
                    )))
                ), dashboardBody((
                  tabsetPanel( id = "tabs",
                               #Information Tab
                               tabPanel(
                                 "Info",
                                 (fluidRow(
                                   includeMarkdown("./introduction.md")
                                 )),
                                 fluidRow(
                                 HTML('<iframe width="560" height="315" ,
                                      src="https://www.youtube.com/embed/th6tD7cBXFM",
                                      frameborder="0" ,
                                      allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" ,
                                      allowfullscreen></iframe>'))
                                 
                               ),
                               #Tutorial Tab
                               tabPanel("Tutorial",
                                        includeMarkdown("./Tutorial/tutorial.md")
                               ),
                               tabPanel("Paper",
                                        htmlOutput("frame")),
                               #Results Tab
                               tabPanel("Methylation Results",
                                        (fluidRow(
                                          box(
                                            width = 12,
                                            title = "Group 3/4 score",
                                            status = "success",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            DTOutput('Mval')
                                          ))),
                                        (fluidRow(
                                          box(
                                            title = "Risk plot",
                                            status = "success",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            plotOutput("figure")),
                                          box(
                                            title = "Selections",
                                            status = "success",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            # h3("Metagene Set:"),
                                            # textOutput("metagenechoice"),
                                            h3("Sample Selected:"),
                                            textOutput("sample"),
                                            h3("Patient's Risk Percentile:"),
                                            textOutput("percentages")
                                          ))),
                                        (fluidRow(
                                          box(
                                            title = "Survival Plot",
                                            status = "success",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            plotOutput("survival")),
                                         box(
                                            title = "Survival Plot",
                                            status = "success",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            plotOutput("survivalage")))),
                                         (fluidRow(
                                           box(
                                             width = 12,
                                             title = "Summary",
                                             status = "success",
                                             solidHeader = TRUE,
                                             collapsible = TRUE,
                                             includeMarkdown("./ResultsSummary.md")
                                           )
                                         ))
                                        ),
                               tabPanel("Expression Results",
                                        (fluidRow(
                                          box(
                                            width = 12,
                                            title = "Group 3/4 score",
                                            status = "warning",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            DTOutput('Mval1')
                                          ))),
                                        (fluidRow(
                                          box(
                                            title = "Risk plot",
                                            status = "warning",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            plotOutput("figureExpression")),
                                          box(
                                            title = "Risk plot",
                                            status = "warning",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                          h3("Sample Selected:"),
                                          textOutput("sampleExpression"),
                                          h3("Patient's Risk Percentile:"),
                                          textOutput("percentagesExpression"),
                                          h3("Patient's survival Percentile:"),
                                          textOutput("survivalExpressionPercentage")
                                          ))),
                                        (fluidRow(
                                          box(
                                            title = "Survival Plot",
                                            status = "warning",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            plotOutput("survivalExpression")),
                                          box(
                                            title = "Survival Plot",
                                            status = "warning",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            plotOutput("survivalageExpression")),
                                        )),
                                        #   box(
                                        #     title = "Risk plot",
                                        #     status = "warning",
                                        #     solidHeader = TRUE,
                                        #     collapsible = TRUE,
                                        #     plotOutput("figure")))),
                                        #   box(
                                        #     title = "Selections",
                                        #     status = "success",
                                        #     solidHeader = TRUE,
                                        #     collapsible = TRUE,
                                        #     h3("Metagene Set:"),
                                        #     textOutput("metagenechoice"),

                                        #   ))),
                                        (fluidRow(
                                          box(
                                            width = 12,
                                            title = "Summary",
                                            status = "warning",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            includeMarkdown("./ResultsSummary.md")
                                          )
                                        ))
                               ),
                               
                               #Download Tab
                               tabPanel("Download",
                                        textInput("filename", "Please insert desired filename", "Group 3/4 Score"),
                                        radioButtons(inputId = "download", label = "Select file type", choices = c("csv"= ".csv", "pdf" = ".pdf")),
                                        downloadButton("down", "Download the results")
                               )
                  )
                )),
                #Font Selection            
                # tags$head(tags$style(HTML('* {font-family: "Courier New"};')))),
  tags$head(tags$style(HTML('* {font-family: "Calibri"};'))))

  
)

server <- function(session, input, output) {
  # Allow larger files to be uploaded
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  #Hide tabs until graph made
  hideTab(inputId = "tabs", target = "Expression Results")
  hideTab(inputId = "tabs", target = "Methylation Results")
  hideTab(inputId = "tabs", target = "Download")
  
  output$frame <- renderUI({
    my_test <- tags$iframe(src="https://www.sciencedirect.com/science/article/pii/S2211124722009718?via%3Dihub", height="900", width="100%")
    my_test
  })
  
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
      
      nmb.mat <- nmb.mat.prepped
 
      # ## interset common genes / probes
      tpms.mat <- match.select(nmb.mat, in.files)
      message(head(tpms.mat))

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

      
      ### define new g3g4 score for projection back onto the original data
      t(rnaseq.H[c(3,1),]) -> g3g4.rnaseq
      message("9")
      apply(g3g4.rnaseq,2,function(x){(1 / (1 + exp(-x)))}) -> logistic.g3g4.rnaseq
      message("10")
      apply(logistic.g3g4.rnaseq,1,function(x){x[2]/(x[1]+x[2])}) -> logistic.g3g4.rnaseq.score
      message("11")
      scaling.function2(logistic.g3g4.rnaseq.score) -> logistic.g3g4.rnaseq.score
      message("12")
      
      ## join and plot the two together (original g3g4 score and g3g4 score projected back onto the same data) kind of a control that it is working
      ## NEED TIDYFT

      message("13")
      
      t(tpms.H[c(3,1),]) -> g3g4.tpms
      message("14")
      apply(g3g4.tpms,2,function(x){(1 / (1 + exp(-x)))}) -> logistic.g3g4.tpms
      message("15")
      apply(logistic.g3g4.tpms,1,function(x){x[2]/(x[1]+x[2])}) -> logistic.g3g4.tpms.score
      message("16")
    
      rbind(logistic.g3g4.rnaseq, logistic.g3g4.tpms) -> scaled.together.logistic.g3g4
      message("list")
      apply( scaled.together.logistic.g3g4,1,function(x){x[2]/(x[1]+x[2])}) ->  scaled.together.logistic.score
      message("apply")
      scaled.together.logistic.score[-1:-length( logistic.g3g4.rnaseq.score)] ->  scaled.together.logistic.score
      message("remove")

      # some times helpful to remove outliers prior to scaling
      outlier.idx <- c(head(order(scaled.together.logistic.score), round((input$outlier)*(length(scaled.together.logistic.score)/100))),
                       tail(order(scaled.together.logistic.score), round((input$outlier)*(length(scaled.together.logistic.score)/100)))
      )
      message("17")
      
      # scale to create final score
      scaling.function3(scaled.together.logistic.score
                       # [-outlier.idx]
                       ) -> logistic.g3g4.tpms.score
      round(logistic.g3g4.tpms.score, digits = 3) -> logistic.g3g4.tpms.score
      as.data.frame(logistic.g3g4.tpms.score) -> logistic.g3g4.tpms.score.df
      message("18")
      
      output$Mval1 <- renderDT (({logistic.g3g4.tpms.score.df
      }),
      options = list(
        pageLength = 10, 
        processing=FALSE),
      selection = "single"
      )
      
      output$figureExpression <-
        renderPlot({
          
          figure.output <-(
            generate_figure_highlight_g3g4Expression(
              logistic.g3g4.tpms.score
              ,input$Mval1_row_last_clicked)
          )
          figure.output
        })
      
      output$percentagesExpression <- renderText (
        
        generate_figure_highlight_g3g4PERC(
          logistic.g3g4.tpms.score,
          input$Mval1_row_last_clicked)
      )
      
      output$sampleExpression <- renderText(input$Mval1_cell_clicked$value) 
      
      output$survivalExpression <- renderPlot({
        figure.output <-(
          survivalcurveplot(
            logistic.g3g4.tpms.score
            ,input$Mval1_row_last_clicked)
        )
        figure.output
      })
      
      output$survivalExpressionPercentage <- renderText({
        figure.output <-(
          survivalcurveplotPERC(
            logistic.g3g4.tpms.score
            ,input$Mval1_row_last_clicked)
        )
        figure.output
      })
      
      output$survivalageExpression <- renderPlot(
        ggplot(df3, aes(x=pred, y=surv, group=age, color = age)) +
          #geom_line() +
          geom_point(alpha = 1/10) +
          geom_line(data = df3.y, aes(x=pred, y=surv, group=age, color = age)) +
          geom_line(data = df3.y, aes(x=pred, y=lo, group=age),linetype="dotted") +
          geom_line(data = df3.y, aes(x=pred, y=up, group=age),linetype="dotted") +
          theme_classic() + xlab("Prediction Metagene") + ylab("Survival") +
          scale_color_manual(values=c('red','dodgerblue')) +
          theme(legend.position = "none") +
          # labs(title = "New plot title", subtitle = "A subtitle") +
          ylim(0,1)
      )

      
      att$done()
      
      showTab(inputId = "tabs", target = "Expression Results", select = TRUE)
      showTab(inputId = "tabs", target = "Download")

      output$down <- downloadHandler(
        filename = function() {
            paste(input$filename,Sys.time(), input$download)
        },
        content ={
          function(file) {
            if (input$download == ".csv")
              write.csv(logistic.g3g4.tpms.score.df,
                        file, row.names = TRUE)
            else
              pdf(file,
                  width = 7,
                  title = paste(input$filename,Sys.time(), ".pdf", sep="_")
              )
            
            grid.table(logistic.g3g4.tpms.score.df)
            
            figure.outputExpression1 <-(
              generate_figure_highlight_g3g4Expression(
                logistic.g3g4.tpms.score
                ,1)
            )
            
            figure.outputSurvival <-(
              survivalcurveplot(
                logistic.g3g4.tpms.score
                ,1)
            )
            
            
            figure.outputAge <- 
              ggplot(df3, aes(x=pred, y=surv, group=age, color = age)) +
              geom_point(alpha = 1/10) +
              geom_line(data = df3.y, aes(x=pred, y=surv, group=age, color = age)) +
              geom_line(data = df3.y, aes(x=pred, y=lo, group=age),linetype="dotted") +
              geom_line(data = df3.y, aes(x=pred, y=up, group=age),linetype="dotted") +
              theme_classic() + xlab("Prediction Metagene") + ylab("Survival") +
              scale_color_manual(values=c('red','dodgerblue')) +
              theme(legend.position = "none") +
              # labs(title = "New plot title", subtitle = "A subtitle") +
              ylim(0,1)
            
            grid.arrange(
              figure.outputExpression1, 
              figure.outputSurvival, 
              figure.outputAge,
              ncol=2)
            
            dev.off()
          }
        })
    })
############################################## METHYLATION ############################################
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
      
      beta2m(temp.processed$betas) -> M.values
      
      ### Round results to 3 figures
      
      metagene <- round(predict(g3.g4.cont.rfe, t(M.values)[,predictors(g3.g4.cont.rfe)]), digits = 3)
      metagene.df <- data.frame('Group.3.4.Score' = metagene)

      output$Mval <- renderDT (
        ({
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
      
      # print(figure.input)

      figure.input <- metagene.df$Group.3.4.Score
      names(figure.input) <- rownames(metagene.df)
      print(figure.input)
      #####
      output$figure <-
        renderPlot({
          
          figure.output <-(
            # figureFile <- "./ecrt20.dist.rds"
            generate_figure_highlight_g3g4(
              figure.input
              ,input$Mval_row_last_clicked)
          )
          figure.output
          })
          # }
      output$percentages <- renderText (
        
        generate_figure_highlight_g3g4PERC(
          figure.input,
          input$Mval_row_last_clicked)
      )
      
      saveRDS(figure.input, file = "./temp/figureinput.RDS")
      
      message(head(df2))
      
      output$survival <- renderPlot({
          figure.output <-(
            survivalcurveplot(
              figure.input
              ,input$Mval_row_last_clicked)
              # c(0.3,0.5),
              # 0.3)
          )
          figure.output
        })
    output$survivalage <- renderPlot(
      ggplot(df3, aes(x=pred, y=surv, group=age, color = age)) +
        #geom_line() +
        geom_point(alpha = 1/10) +
        geom_line(data = df3.y, aes(x=pred, y=surv, group=age, color = age)) +
        geom_line(data = df3.y, aes(x=pred, y=lo, group=age),linetype="dotted") +
        geom_line(data = df3.y, aes(x=pred, y=up, group=age),linetype="dotted") +
        # theme_classic() + 
        xlab("Prediction Metagene") + 
        ylab("Survival") +
        scale_color_manual(values=c('red','dodgerblue')) +
        # labs(title = "New plot title", subtitle = "A subtitle") +
        ylim(0,1)
    )


          # https://rstudio.github.io/DT/shiny.html


      
    #############
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
      
      showTab(inputId = "tabs", target = "Methylation Results", select = TRUE)
      showTab(inputId = "tabs", target = "Download")
      
      #For displaying currently selected sample
      
      output$sample <- renderText(input$Mval_cell_clicked$value) 
      # figure.input <- readRDS("./temp/figureinput.RDS")
      
      figinput <- readRDS("./temp/figureinput.RDS")
      figureinputDF <- as.data.frame(figinput)


      output$down <- downloadHandler(
        filename = function() {
          paste(input$metagenes,Sys.time(), input$download, sep=".")
        },
        content ={ 
          function(file) {
            if (input$download == "csv")
              
            
              write.csv(figureinputDF,
                        file, row.names = TRUE)

            ####### if not csv then write a pdf
            else
              pdf(file,
                  width = 14,
                  title = paste("output")
                  )
            
            colnames(figureinputDF) <- "G3/G4_Score"
            grid.table(figureinputDF)


            figure.output <-(
              generate_figure_highlight_g3g4(
                figure.input
                ,input$Mval_row_last_clicked)
            )

            figure.output1 <-(
              survivalcurveplot(
                figure.input
                ,input$Mval_row_last_clicked)
            )

            
            figure.output2 <- 
              ggplot(df3, aes(x=pred, y=surv, group=age, color = age)) +
              #geom_line() +
              geom_point(alpha = 1/10) +
              geom_line(data = df3.y, aes(x=pred, y=surv, group=age, color = age)) +
              geom_line(data = df3.y, aes(x=pred, y=lo, group=age),linetype="dotted") +
              geom_line(data = df3.y, aes(x=pred, y=up, group=age),linetype="dotted") +
              # theme_classic() + 
              xlab("Prediction Metagene") + 
              ylab("Survival") +
              scale_color_manual(values=c('red','dodgerblue')) +
              # labs(title = "New plot title", subtitle = "A subtitle") +
              ylim(0,1)
            
            # plot(figure.output2)
            
            grid.arrange(
              figure.output, 
              figure.output1, 
              figure.output2,
              ncol=2)

              dev.off()
            
          }

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

}

shinyApp(server = server, ui = ui)