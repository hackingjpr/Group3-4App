source("./AppSourceFunctions1.10.R")

message("Begin the app")
ui <- shiny::fluidPage(
  message("UI initiating"),
  #Loading bar
  useAttendant(),
  #Dashboard looks etc
  dashboardPage(
    dashboardHeader(title = "", dropdownMenuOutput("messageMenu")),
    dashboardSidebar(
      #logo
      # (img(src='Free_Sample_By_Wix (5).jpg', align = "center")),
      (img(src = 'Free_Sample_By_Wix%20(12).jpg', align = "center")),
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
          accept = c(".rds", ".csv", ".txt")
        )
      ),
      conditionalPanel(
        condition = "input.expmeth == 'exp'",
        fileInput(
          "expFile",
          "Upload Expression files:",
          multiple = TRUE,
          accept = c(".rds", ".csv", ".txt")
        )
      ),
      conditionalPanel(
        condition = "input.expmeth == 'exp'",
        radioButtons(
          "scaling",
          "Scale using Williamson Et Al dataset or use your uploaded data?",
          c(
            "Williamson et al Dataset" = "ours",
            "Uploaded Dataset" = "yours"
          )
        )
      ),
      conditionalPanel(
        condition = "input.scaling == 'yours'",
        sliderInput(
          "outlier",
          "OUTLIER SELECTION WHAT SHOULD I CALL THIS??",
          min = 1,
          max = 4,
          value = 1,
          step = 1,
          round = 0,
          animate = FALSE
        )
      )
      ,
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
          color = "warning",
          style = "unite"
        )
      ),
      conditionalPanel(
        condition = "input.expmeth == 'meth'",
        actionBttn(
          inputId = "bttn2",
          label = "Generate Group 3/4 Score, Methylation",
          color = "success",
          style = "unite"
        )
      ),
      
      #Reset button
      column(width = 12,
             "Want to upload a different data set?"),
      column(
        width = 12,
        actionBttn(
          inputId = "bttn3",
          label = "Reset",
          color = "danger",
          style = "unite",
          size = "lg"
        ),
        align = "center"
      ),
      #Loading Bar
      attendantBar("progress-bar"),
      hr(),
      tags$footer(HTML(
        paste(
          "Author: James Hacking,",
          "<br/>",
          "Date Created: 26-07-2022,",
          "<br/>",
          "Copyright (c) James Hacking, 2022,",
          "<br/>",
          h5(
            a("Email: james.hacking@ncl.ac.uk", href = "mailto:james.hacking@ncl.ac.uk")
          ),
          h4(
            "Disclaimer : This app is designed exclusively for research purposes and is strictly not for diagnostic use."
          ),
          h5("Supported by funding from LoveOliver, Children with Cancer UK and CRUK.")
          # shinythemes::themeSelector,()
          
        )
      ))
    ),
    dashboardBody((
      tabsetPanel(
        id = "tabs",
        #Information Tab
        tabPanel(style = "padding-left:15px",
                 "Info",
                 (fluidRow(
                   includeMarkdown("./introduction.md")
                 )),
                 fluidRow(
                   HTML(
                     '<iframe width="560" height="315" ,
                                      src="https://www.youtube.com/embed/th6tD7cBXFM",
                                      frameborder="0" ,
                                      allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" ,
                                      allowfullscreen></iframe>'
                   )
                 )),
        #Tutorial Tab
        tabPanel("Tutorial",
                 includeMarkdown("./Tutorial/tutorialIO.md")),
        tabPanel("Paper",
                 htmlOutput("frame")),
        #Results Tab
        tabPanel(
          "Methylation Results",
          (fluidRow(
            box(
              width = 12,
              title = "Group 3/4 score",
              status = "success",
              solidHeader = TRUE,
              collapsible = TRUE,
              DTOutput('Mval')
            )
          )),
          (fluidRow(
            box(
              title = "Group 3/4 Plot",
              status = "success",
              solidHeader = TRUE,
              collapsible = TRUE,
              plotOutput("figure")
            ),
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
              textOutput("percentages"),
              
              h3("Patient's Survival Percentile:"),
              textOutput("survivalPercentage"),
              
              h4("Disclaimer : This app is designed exclusively for research purposes and is strictly not for diagnostic use.")
            )
          )),
          (fluidRow(
            box(
              title = "Survival Plot: No Risk Factors Considered",
              status = "success",
              solidHeader = TRUE,
              collapsible = TRUE,
              plotOutput("survival")
            ),
            box(
              title = "Survival Plot: Age Considered",
              status = "success",
              solidHeader = TRUE,
              collapsible = TRUE,
              plotOutput("survivalage")
            )
          )),
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
        tabPanel(
          "Expression Results",
          (fluidRow(
            box(
              width = 9,
              title = "Group 3/4 score",
              status = "warning",
              solidHeader = TRUE,
              collapsible = TRUE,
              DTOutput('Mval1')
            )
          )),
          (fluidRow(
            box(
              title = "Group 3/4 Plot",
              status = "warning",
              solidHeader = TRUE,
              collapsible = TRUE,
              plotOutput("figureExpression")
            ),
            box(
              title = "Selected Sample Information",
              width = 3,
              status = "warning",
              solidHeader = TRUE,
              collapsible = TRUE,
              h3("Sample Selected:"),
              textOutput("sampleExpression"),
              # h3("Patients Group 3/4 Score"),
              # textOutput("scoreExpression"),
              h3("Patient's Risk Percentile:"),
              textOutput("percentagesExpression"),
              h3("Patient's Survival Percentile:"),
              textOutput("survivalExpressionPercentage"),
              
              tags$b(
                "Disclaimer : This app is designed exclusively for research purposes and is strictly not for diagnostic use.")
            )
          )),
          (fluidRow(
            box(
              title = "Survival Plot: No Risk Factors Considered",
              status = "warning",
              solidHeader = TRUE,
              collapsible = TRUE,
              plotOutput("survivalExpression")
            ),
            box(
              title = "Survival Plot: Age Considered",
              status = "warning",
              solidHeader = TRUE,
              collapsible = TRUE,
              plotOutput("survivalageExpression")
            ),
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
        tabPanel(
          "Download Expression",
          textInput("filename", "Please insert desired filename", "Group 3/4 Score"),
          radioButtons(
            inputId = "downloadExp",
            label = "Select file type",
            choices = c("csv" = ".csv", "pdf" = ".pdf")
          ),
          downloadButton("downExpression", "Download the results")
        ),
        tabPanel(
          "Download Methylation",
          textInput("filename", "Please insert desired filename", "Group 3/4 Score"),
          radioButtons(
            inputId = "downloadMethy",
            label = "Select file type",
            choices = c("csv" = ".csv", "pdf" = ".pdf")
          ),
          downloadButton("downMethylation", "Download the results")
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
  hideTab(inputId = "tabs", target = "Download Expression")
  hideTab(inputId = "tabs", target = "Download Methylation")
  
  
  
  output$frame <- renderUI({
    my_test <-
      tags$iframe(src = "https://www.sciencedirect.com/science/article/pii/S2211124722009718?via%3Dihub",
                  height = "900",
                  width = "100%")
    my_test
  })
  
  #Loading bar
  att <- Attendant$new("progress-bar")
  unlink("./temp/", recursive = T)
  if (!dir.exists("./temp")) {
    dir.create("./temp")
  }
  #fix the overwriting issue
  fs::path("./temp/", createRandString()) -> tempDIR
  dir.create(tempDIR)
  ############################################## EXPRESSION ############################################## 
  Mvals <-
    observeEvent(
      input$bttn1,
      {
        req(input$expFile)
        att$set(0, text = "Loading") #Start at 10%
        att$auto(ms = 1600, value = 0.0048) # automatically increment
        message("Loading Bar")
        
        
        
        require(R.utils)
        
        input$expFile$datapath -> in.files
        message(in.files)
        paste0(tempDIR, "/", input$expFile$name) -> out.files
        message(out.files)
        message("TempDir")
        
        copied <- file.copy(in.files, out.files)
        message("copy")
        

        input.file <- in.files
        
        if (file_ext(input.file) == "rds") {
          in.files <- readRDS(file = input.file)
        } else if (file_ext(input.file) == "csv") {
          in.files <- read.csv(file = input.file, row.names = 1)
        } else if (file_ext(input.file) == "txt") {
          in.files <- read.delim(file = input.file)
        } else {
          message("file not right format!")
        }
        
        
        message("File loaded")
        nmb.mat <- nmb.mat.prepped
        
        # saveRDS(in.files, "~/Group3-4App/temp/csvfile.rds")
        # in.files <- readRDS("~/Group3-4App/temp/csvfile.rds")
        
        # ## interset common genes / probes
        tpms.mat <- match.select(nmb.mat, in.files)
        # message(head(tpms.mat))
        #saveRDS(tpms.mat, "~/Group3-4App/temp/tpms.mat10.rds")
        
        
        message("6")
        
        
        
        ## project using pseudo-inverse & post-projection normalise
        # project back onto the same dataset
        rnaseq.H <- project.NMF(input.array = nmb.mat,
                                nmf.result = nmf.res)
        message("7")
        
        tpms.matrix <- as.matrix(tpms.mat)
        if (ncol(tpms.matrix) == 1) {
          colnames(in.files) -> colnames(tpms.matrix)
        }
        message(object.size(tpms.matrix))
        message(dim(tpms.matrix))
        # project onto fresh dataset
        tpms.H <- project.NMF(input.array = tpms.matrix,
                              nmf.result = nmf.res)
        message("8")
        
        
        ### define new g3g4 score for projection back onto the original data
        t(rnaseq.H[c(3, 1), ]) -> g3g4.rnaseq
        message("9")
        apply(g3g4.rnaseq, 2, function(x) {
          (1 / (1 + exp(-x)))
        }) -> logistic.g3g4.rnaseq
        message("10")
        apply(logistic.g3g4.rnaseq, 1, function(x) {
          x[2] / (x[1] + x[2])
        }) -> logistic.g3g4.rnaseq.score
        message("11")
        if (input$scaling == "ours")
          scaling.function3(logistic.g3g4.rnaseq.score) -> logistic.g3g4.rnaseq.score
        
        else
          scaling.function(logistic.g3g4.rnaseq.score) -> logistic.g3g4.rnaseq.score
        message("12")
        
        ## join and plot the two together (original g3g4 score and g3g4 score projected back onto the same data) kind of a control that it is working
        ## NEED TIDYFT
        
        message("13")
        
        t(tpms.H[c(3,1),]) -> g3g4.tpms
        
        message("14")
        
        apply(g3g4.tpms, 2, function(x) {
          (1 / (1 + exp(-x)))
        }) -> logistic.g3g4.tpms
        
        message("15")
        
        if(is.null(dim(logistic.g3g4.tpms))){
          logistic.g3g4.tpms[2] / (logistic.g3g4.tpms[1] + logistic.g3g4.tpms[2]) -> logistic.g3g4.tpms.score
          apply(logistic.g3g4.tpms, 1, function(x) {
            x[2] / (x[1] + x[2])
          }) ->  logistic.g3g4.tpms.score
          message("is.null(dim(logistic.g3g4.tpms))")
          scaling.function3(logistic.g3g4.tpms.score) -> logistic.g3g4.tpms.score
        }else{
          apply(logistic.g3g4.tpms, 1, function(x) {
            x[2] / (x[1] + x[2])
          }) -> logistic.g3g4.tpms.score
          
          # evaluate this expression if(input$outlier !=0 & round((5)*(length(scaled.together.logistic.score)/100))){ xxxxxxx    }
          # some times helpful to remove outliers prior to scaling
          
          # some times helpful to remove outliers prior to scaling
          # outlier.idx <- c(head(order(scaled.together.logistic.score), round((input$outlier)*(length(scaled.together.logistic.score)/100))),
          #                  tail(order(scaled.together.logistic.score), round((input$outlier)*(length(scaled.together.logistic.score)/100)))
          # )
          
          mean(logistic.g3g4.tpms.score) -> mean.logistic.g3g4.tpms.score
          message(mean.logistic.g3g4.tpms.score)
          sd(logistic.g3g4.tpms.score) -> sd.logistic.g3g4.tpms.score
          message(sd.logistic.g3g4.tpms.score)
          
          if (input$outlier == 0){
            upper.limit <- 1
            lower.limit <- 0
          }
          else{
            upper.limit <-
              ((input$outlier) * sd.logistic.g3g4.tpms.score) +  mean.logistic.g3g4.tpms.score
            lower.limit <-
              mean.logistic.g3g4.tpms.score - ((input$outlier) * sd.logistic.g3g4.tpms.score)
          }
          
          message(upper.limit)
          message(lower.limit)
          message(input$outlier)
          
          outlier.idx <-
            which(logistic.g3g4.tpms.score > upper.limit |
                    logistic.g3g4.tpms.score < lower.limit)
          
          message(length(outlier.idx))
          message(outlier.idx)
          
          if (length(outlier.idx) != 0 & input$scaling == "yours") {
            apply(logistic.g3g4.tpms, 1, function(x) {
              x[2] / (x[1] + x[2])
            }) ->  logistic.g3g4.tpms.score
            message("length(outlier.idx) != 0 & input$scaling == 'yours'")
            removed <- logistic.g3g4.tpms.score[-outlier.idx]
            message(removed)
            # scaling.function(logistic.g3g4.tpms.score
            #                  [-outlier.idx]) -> logistic.g3g4.tpms.score
            scaling.function(removed) -> logistic.g3g4.tpms.score
            message(logistic.g3g4.tpms.score)
          }else if (length(outlier.idx) == 0 & input$scaling == "yours") {
            apply(logistic.g3g4.tpms, 1, function(x) {
              x[2] / (x[1] + x[2])
            }) ->  logistic.g3g4.tpms.score
            message("length(outlier.idx) == 0 & input$scaling == 'yours'")
            scaling.function(logistic.g3g4.tpms.score) -> logistic.g3g4.tpms.score
          } else{
            apply(logistic.g3g4.tpms, 1, function(x) {
              x[2] / (x[1] + x[2])
            }) ->  logistic.g3g4.tpms.score
            scaling.function3(logistic.g3g4.tpms.score) -> logistic.g3g4.tpms.score
          }
        }
        
        round(logistic.g3g4.tpms.score, digits = 3) -> logistic.g3g4.tpms.score
        data.frame('Group.3.4.Score' = logistic.g3g4.tpms.score) -> logistic.g3g4.tpms.score.df
        # logistic.g3g4.tpms.score.df <- ('Group.3.4.Score' = logistic.g3g4.tpms.score)
        message("18")
        
        output$Mval1 <- renderDT (({
          logistic.g3g4.tpms.score.df
        }),
        options = list(pageLength = 10,
                       processing = FALSE),
        selection = "single")
        
        message(is.null(input$Mval1_row_last_clicked))
        message("input$Mval1_row_last_clicked")
        
        # saveRDS(logistic.g3g4.tpms.score, "./temp/logistic.rds")
        # saveRDS(input$Mval1_row_last_clicked, "./temp/ROW.RDS")
        output$figureExpression <-
          renderPlot({
            figure.output <- (
              generate_figure_highlight_g3g4Expression(logistic.g3g4.tpms.score
                                                       , input$Mval1_row_last_clicked)
            )
            figure.output
          })
        
        # generate_figure_highlight_g3g4Expression(logistic, ROW)
        
        output$percentagesExpression <- renderText ({
          paste(
            generate_figure_highlight_g3g4PERC(logistic.g3g4.tpms.score,
                                               input$Mval1_row_last_clicked),
            "of patients had a lower group 3/4 score than this."
          )
        })
        
        output$sampleExpression <- renderText(
          if (is.null(input$Mval1_row_last_clicked))  {
            row.names(logistic.g3g4.tpms.score.df) [1]
          } 
          else {
            input$Mval1_cell_clicked$value
          }
        )
        
        # input$Mval1_cell_clicked$value) 
        
        output$survivalExpression <- renderPlot({
          figure.output <-(
            survivalcurveplot(
              logistic.g3g4.tpms.score
              ,input$Mval1_row_last_clicked)
          )
          figure.output
        })
        
        # survivalcurveplot(logistic,ROW)
        
        output$survivalExpressionPercentage <- renderText({
          figure.output <-(
            survivalcurveplotPERC(
              logistic.g3g4.tpms.score
              ,input$Mval1_row_last_clicked)
          )
          figure.outputOLDER <-(
            SurvivalAgePlotPerc5(
              logistic.g3g4.tpms.score
              ,input$Mval1_row_last_clicked)
          )
          figure.outputYOUNGER <-(
            SurvivalAgePlotPerc4(
              logistic.g3g4.tpms.score
              ,input$Mval1_row_last_clicked)
          )
          
          paste("According to the dataset in Williamson et als retrospective survival cohort, patients with this score on average had a 5 year survival percentage of",
                figure.output,
                ". This does not take into account other risk factors. If age is taken into account the percentage would change to",
                figure.outputYOUNGER, "if the patient was under three years old and",
                figure.outputOLDER, "if older than three. (See below)")
        })
        
        output$survivalageExpression <- renderPlot(
          SurvivalAgePlot(logistic.g3g4.tpms.score,
                          input$Mval1_row_last_clicked)
        )
        
        
        att$done()
        
        showTab(inputId = "tabs", target = "Expression Results", select = TRUE)
        showTab(inputId = "tabs", target = "Download Expression")
        
        output$downExpression <- downloadHandler(
          filename = function() {
            paste(input$filename,Sys.time(), input$downloadExp)
          },
          content ={
            function(file) {
              if (input$downloadExp == ".csv")
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
                SurvivalAgePlot(logistic.g3g4.tpms.score,
                                1)
              
              grid.arrange(
                figure.outputExpression1, 
                figure.outputSurvival, 
                figure.outputAge,
                ncol=2)
              
              dev.off()
            }
          })
      })
  ############################################## METHYLATION ############################################## 
  observeEvent(input$bttn2, {
    att$set(10, text = "Loading") #Start at 10% 
    att$auto(ms = 1600, value = 0.01) # automatically increment
    
    #need idat file to run
    req(input$methFile)
    require(R.utils)
    

    
    input$methFile$datapath -> in.files
    message(in.files)
    paste0(tempDIR, "/", input$methFile$name) -> out.files
    message(out.files)
    message("TempDir")
    
    copied <- file.copy(in.files, out.files)
    message("copy")
    
    
    input.file <- in.files
    
    if (file_ext(input.file) == "rds") {
      in.files <- readRDS(file = input.file)
    } else if (file_ext(input.file) == "csv") {
      in.files <- read.csv(file = input.file, row.names = 1)
    } else if (file_ext(input.file) == "txt") {
      in.files <- read.delim(file = input.file, row.names = 1)
    } else {
      message("file not right format!")
    }
    
    # temp.base <- get_basenames(tempDIR)
    
    #saveRDS(temp.base, "./temp/temp.base.rds")
    
    
    cat("Timing start\n")
    ptm <- proc.time()
    
    # temp.processed <- process_idats(input.file)
    # on.exit({
    #   att$done()
    # })
    
    # saveRDS(temp.processed$betas, "./Inputs/temp-processed-betas.rds")
    # write.csv(temp.processed$betas, "./Inputs/temp.processed.betas.csv")
    # temp.processeddf <- as.data.frame(temp.processed$betas)
    # write.csv(temp.processeddf, "./Inputs/temp.processed.betas.txt")

    
    # READING IN M.VALS
    message("I got to beta2m")
    beta2m(in.files) -> M.values
    
    if(ncol(M.values)==1){
      ### for single sample
      t(data.frame(t(M.values)[,predictors(g3.g4.cont.rfe)])) -> input.df
      colnames(M.values) -> rownames(input.df)
    }else{
      t(M.values)[,predictors(g3.g4.cont.rfe)] -> input.df
    }
    
    ### Round results to 3 figures
    metagene <- round(predict(g3.g4.cont.rfe, input.df), digits = 3)
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
    
    
    figure.input <- metagene.df$Group.3.4.Score
    names(figure.input) <- rownames(metagene.df)
    print(figure.input)
    #####
    output$figure <-
      renderPlot({
        
        figure.output <-(
          generate_figure_highlight_g3g4(
            figure.input
            ,input$Mval_row_last_clicked)
        )
        figure.output
      })
    # }
    output$percentages <- renderText (
      
      {
        paste(
          generate_figure_highlight_g3g4PERC(
            figure.input,
            input$Mval_row_last_clicked),
          "of patients had a lower group 3/4 score than this.")
      }
    )
    
    saveRDS(figure.input, file = "./temp/figureinput.RDS")
    
    
    message(head(df2))
    
    
    
    output$survival <- renderPlot({
      figure.output <-(
        survivalcurveplot(
          figure.input
          ,input$Mval_row_last_clicked)
      )
      figure.output
    })
    
    output$survivalPercentage <- renderText({
      figure.output <-(
        survivalcurveplotPERC(
          figure.input
          ,input$Mval_row_last_clicked))
      figure.outputOLDER <-(
        SurvivalAgePlotPerc5(
          figure.input
          ,input$Mval1_row_last_clicked)
      )
      figure.outputYOUNGER <-(
        SurvivalAgePlotPerc4(
          figure.input
          ,input$Mval1_row_last_clicked)
      )
      
      paste("According to the dataset in Williamson et als retrospective survival cohort, patients with this score on average had a 5 year survival percentage of",
            figure.output,
            ". This does not take into account other risk factors. If age is taken into account the percentage would change to",
            figure.outputYOUNGER, "if the patient was under three years old and",
            figure.outputOLDER, "if older than three. (See below)")
      
    })
    
    
    output$survivalage <- renderPlot(
      SurvivalAgePlot(figure.input,
                      input$Mval_row_last_clicked)
    )
    
    # https://rstudio.github.io/DT/shiny.html
    
    #####
    # output$metagenechoice <- renderText({input$metagenes})
    
    showTab(inputId = "tabs", target = "Methylation Results", select = TRUE)
    showTab(inputId = "tabs", target = "Download Methylation")
    
    #For displaying currently selected sample
    
    output$sample <- renderText(input$Mval_cell_clicked$value) 
    # figure.input <- readRDS("./temp/figureinput.RDS")
    
    figinput <- readRDS("./temp/figureinput.RDS")
    figureinputDF <- as.data.frame(figinput)
    colnames(figureinputDF) <- "G3/G4 Score"
    
    
    
    output$downMethylation <- downloadHandler(
      filename = function() {
        paste(input$filename,Sys.time(), input$downloadMethy, sep=".")
      },
      content ={ 
        function(file) {
          if (input$downloadMethy == ".csv")
            # write.csv(figureinputDF,
            #           file, row.names = TRUE)
            write.csv(figureinputDF,
                      file, row.names = TRUE)
          # read.csv("./temp/figureinputDF.csv")
          
          ####### if not csv then write a pdf
          else
            pdf(file,
                width = 14,
                title = paste(input$filename,Sys.time(), input$downloadMethy, sep=".")
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
            SurvivalAgePlot(figure.input,
                            1)
          
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