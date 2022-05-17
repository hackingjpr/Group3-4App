source("./source_functions.R")


ui <- shiny::fluidPage(
  
  # tags$h1("NAME OF APP TO FOLLOW"),
  #Loading bar
  useAttendant(),
  #Dashboard looks etc
  dashboardPage(dashboardHeader(title= "", dropdownMenuOutput("messageMenu")), 
                dashboardSidebar(
                  #logo
                  # (img(src='Free_Sample_By_Wix (5).jpg', align = "center")),
                  (img(src='Free_Sample_By_Wix%20(5).jpg', align = "center")),
                  #uploading idat
                  h2("Step 1"),
                  fileInput(
                    "idatFile",
                    "Upload idat file:",
                    multiple = TRUE,
                    accept = ".idat"),
                  #Step 2
                  h2("Step 2"),
                  selectInput("metagenes", "Select Metagene Set", c("MRT (ATRT & ECRT)", "ATRT", "ECRT")),
                  #set it away
                  h2("Step 3"),
                  actionBttn(
                    inputId = "bttn1",
                    label = "Generate Risk Values",
                    color = "primary",
                    style = "fill"
                  ),
                  
                  #Reset button
                  actionBttn(
                    inputId = "bttn2",
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
                      "Date Created: 28-04-2022,", "<br/>", 
                      "Copyright (c) James Hacking, 2022,", "<br/>",
                      h5(a("Email: james.hacking@ncl.ac.uk", href="mailto:james.hacking@ncl.ac.uk"))
                      
                    )))
                ), dashboardBody((
                  tabsetPanel( id = "tabs",
                               #Information Tab
                               tabPanel(
                                 "Info",
                                 (fluidRow(
                                   
                                   "Here is where I will talk about what you can do and what it is all about etc.",
                                   h1("THIS IS FOR RESEARCH PURPOSES ONLY, DO NOT USE FOR DIAGNOSTICS")
                                   
                                 )
                                 )
                               ),
                               #Tutorial Tab
                               tabPanel("Tutorial",
                                        fluidRow("Tutorial here")),
                               #Results Tab
                               tabPanel("Results Table",
                                        (fluidRow(
                                          # column(100,
                                          box(
                                            width = 12,
                                            title = "Risk Values",
                                            status = "info",
                                            solidHeader = TRUE,
                                            DTOutput('Mval')
                                            # )
                                          ))),
                                        # (fluidRow(plotOutput("figure")%>% withSpinner(color="#0dc5c1"))),
                                        (fluidRow(
                                          box(
                                            title = "Risk plot", 
                                            status = "warning", 
                                            solidHeader = TRUE,
                                            collapsible = FALSE,
                                            plotOutput("figure")),
                                          box(
                                            title = "Selections",
                                            status = "success",
                                            solidHeader = TRUE,
                                            "Metagene Set:",
                                            textOutput("metagenechoice"),
                                            "Sample Selected:",
                                            textOutput("sample")
                                          ))
                                         #textOutput("time")
                                        )),
                               
                               #Download Tab
                               tabPanel("Download",
                                        textInput("filename", "Please insert desired filename", "M-values"),
                                        radioButtons(inputId = "download", label = "Select file type", choices = c("csv", "pdf")),
                                        downloadButton("down", "Download the results"),
                               )
                  )
                )),
                #Font Selection            
                tags$head(tags$style(HTML('* {font-family: "Courier New"};')))),
  
  # setBackgroundColor(
  #   #color = c("#F7FBFF", "#1E90FF"),
  #   color = c("#FFFFFF", "#FFFFFF"),
  #   gradient = "radial",
  #   direction = c("top", "left")
  # ),
)

server <- function(session, input, output) {
  # Allow larger files to be uploaded
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  #Hide tabs until graph made
  hideTab(inputId = "tabs", target = "Results Table")
  hideTab(inputId = "tabs", target = "Download")
  
  #Loading bar
  att <- Attendant$new("progress-bar")
  
  if(!dir.exists("./temp")){
    dir.create("./temp")
  }
  #fix the overwriting issue
  fs::path("./temp/",createRandString()) -> tempDIR
  dir.create(tempDIR)
  
  Mvals <-
    observeEvent(input$bttn1, {
      att$set(10, text = "Loading") #Start at 10% 
      att$auto() # automatically increment
      
      #need idat file to run
      req(input$idatFile)
      require(R.utils)
      
      
      input$idatFile$datapath -> in.files
      message(in.files)
      paste0(tempDIR,"/", input$idatFile$name) -> out.files
      message(out.files)
      
      copied <- file.copy(in.files, out.files)
      
      
      temp.base <- get_basenames(tempDIR)
      
      cat("Timing start\n")
      ptm <- proc.time()
      
      temp.processed <- process_idats(temp.base)
      on.exit({
        att$done()
      })
      
      metagene.react <-  
        
        if (input$metagenes == "MRT (ATRT & ECRT)") {
          ALL -> meta
        }
      if (input$metagenes == "ATRT") {
        ATRT -> meta
      }
      if (input$metagenes == "ECRT") {
        ECRT -> meta
      }
      #this needs to be reactive
      # output$test.res <- reactive({extract.metagene(
      test.res <- extract.metagene(
        as.character(meta[[1]]$genes),
        as.numeric(meta[[1]]$weights),
        beta2m(temp.processed$betas),
        as.numeric(meta[[2]])
      )
    # )})
      round(test.res, digits = 3) -> test.res
      
      output$Mval <- renderDT (({test.res
      }),
      options = list(
        pageLength = 10,
        # initComplete = I("function(settings, json) {alert('Done.');}"),
        processing=FALSE),
      #selection = list(target = 'row+column'),
      selection = "single"
      # target = "row+column",
      #selection = list(selected = c(1))
      )
      
      print(test.res)
      figure.input <- test.res$Risk_Value
      names(figure.input) <- rownames(test.res)
      print(figure.input)
      
      att$done(text = "Complete")
      
      #indexRow <- reactive({input$Mval_row_last_clicked})
      #print(input$Mval_row_last_clicked)
      
      
      
      output$time <- renderText({proc.time() - ptm})
      # cat(test)
      # test <- input$Mval_rows_all
      # test <- input$Mval_cell_clicked
      # s1 <- test$col
      # s2 <- test$value
      rowSelect <- reactive({input$Mval_rows_selected})
      
      output$figure <-
        renderPlot(
          # test <- input$Mval_rows_all,
          
          # https://rstudio.github.io/DT/shiny.html
          
          
          
          generate_figure_highlight(
            # c(#input$Mval_cells_selected,
            #  #"test$col" = test,
            # #input$Mval_cell_clicked$value,
            # "1" = 1,
            # "0.4" = 0.4)
            figure.input
            ,input$Mval_row_last_clicked)
          # (c("patientA" = 0.1,
          #    test,
          #    "fljkfd" = 1))
        )
      output$metagenechoice <- renderText({input$metagenes})
      showTab(inputId = "tabs", target = "Results Table", select = TRUE)
      showTab(inputId = "tabs", target = "Download")
      
      #For displaying currently selected sample
      # output$sample <- renderText({
      # s = input$Mval_cell_clicked$value
      # if (length(s)) {
      #   cat('These Samples were selected:\n\n')
      #   cat(s, sep = ', ')
      # }})
      
      output$sample <- renderText(input$Mval_cell_clicked$value) 
      #output$sample <- renderText(input$Mval_row_last_clicked) 
      
      #Message popup
      # from <- "Server"
      # message <- "Processing complete!"
      # messageData <- data.frame(from, message)
      # output$messageMenu <- renderMenu({
      #   # Code to generate each of the messageItems here, in a list. This assumes
      #   # that messageData is a data frame with two columns, 'from' and 'message'.
      #   msgs <- apply(messageData, 1, function(row) {
      #     messageItem(from = row[["from"]], message = row[["message"]])
      #   })
      #   
      #   # This is equivalent to calling:
      #   #   dropdownMenu(type="messages", msgs[[1]], msgs[[2]], ...)
      #   dropdownMenu(type = "messages", badgeStatus = "success", .list = msgs)
      # })
    })
  # output$figure <- renderPlot(generate_figure
  #                             #(input$input$tableId_row_last_clicked)
  #                             (c("patientA" = 0.1,
  #                                input$input$tableId_row_last_clicked,
  #                                "fljkfd" = 1))
  # )
  #Reset session and delete tempDIR
  observeEvent(input$bttn2, {
    unlink(tempDIR, recursive = T)
    session$reload()
    return()
    print("session reload not working")
  })
  renderText(output$metagenes <- input$metagenes)
  output$down <- downloadHandler(
    filename = function() {
      paste(input$filename,Sys.time(), ".csv", sep="_")
    },
    content ={ 
      function(file) {
        write.csv(Mvals(), file, row.names = TRUE)
      }
      # if (input$download == "csv")
      #   write.csv(output$Mval)
      # else
      #   pdf(output$Mval,
      #       width = 14
      #   )
      # dev.off()
    })
  
  
  
  
  
  
}

shinyApp(server = server, ui = ui)
