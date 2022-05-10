source("./source_functions.R")


ui <- shiny::fluidPage(
  # tags$h1("NAME OF APP TO FOLLOW"),
  # 
  dashboardPage(dashboardHeader(title= "Name of App", dropdownMenuOutput("messageMenu")), 
                dashboardSidebar(
                  (img(src='Free_Sample_By_Wix.jpg', align = "center")),
                  fileInput(
    "idatFile",
    "Upload idat file:",
    multiple = TRUE,
    accept = ".idat"),
    selectInput("metagenes", "Metagenes", c("ALL", "ATRT", "ECRT")),
    # column(
    #   3,
    actionBttn(
      inputId = "bttn1",
      label = "Generate M Values",
      color = "primary",
      style = "fill"
    ),
    
    # fluidRow( 
    actionBttn(
      inputId = "bttn2",
      label = "Reset",
      color = "danger",
      style = "fill",
      size = "sm"
    ),
    hr(),
    tags$footer( HTML(
      paste(
      "Author: James Hacking,", "<br/>",  
                             "Date Created: 28-04-2022,", "<br/>", 
                             "Copyright (c) James Hacking, 2022,", "<br/>",
              "Email: james.hacking@ncl.ac.uk"
    )))
    ), dashboardBody((
      tabsetPanel(
        tabPanel(
          "Info",
          (fluidRow(
            
            "Here is where I will talk about what you can do and what it is all about etc."

          )
          )
        ),
        tabPanel("Results Table",
                 (fluidRow(
                   textOutput("metagenes"),
                   column(12,
                          DTOutput('Mval')))),
                  (fluidRow(plotOutput("figure"))),
                  textOutput("time")
                 ),
        
        
        tabPanel("Download",
                 textInput("filename", "Please insert desired filename", "M-values"),
                 radioButtons(inputId = "download", label = "Select file type", choices = c("csv", "pdf")),
                 downloadButton("down", "Download the results"),
        )
      )
    )),
                tags$head(tags$style(HTML('* {font-family: "Courier New"};')))),
  
  setBackgroundColor(
    #color = c("#F7FBFF", "#1E90FF"),
    color = c("#FFFFFF", "#FFFFFF"),
    gradient = "radial",
    direction = c("top", "left")
  ),
    )
  
server <- function(session, input, output) {
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  Mvals <-
    observeEvent(input$bttn1, {
      req(input$idatFile)
      require(R.utils)
      fs::path(tempdir(),createRandString()) -> tempDIR
      dir.create(tempDIR)
      
      withProgress(message = 'beep boop, doing basey things', value = 1, {
        input$idatFile$datapath -> in.files
        message(in.files)
        paste0(tempDIR,"/", input$idatFile$name) -> out.files
        message(out.files)

        copied <- file.copy(in.files, out.files)

        
        temp.base <- get_basenames(tempDIR)})
      withProgress(message = "beep boop, doing idat things, trust me it will finish I just haven't worked out how to link the loading bar to progress", value = 0, {
        
        cat("Timing start\n")
        ptm <- proc.time()
        # for (i in 1:150){
        # incProgress(amount = 1/150)
        #   Sys.sleep(1)}
        temp.processed <- process_idats(temp.base)
      })
      output$Mval <- renderDT (({
        if (input$metagenes == "ALL") {
          ALL -> meta
        }
        if (input$metagenes == "ATRT") {
          ATRT -> meta
        }
        if (input$metagenes == "ECRT") {
          ECRT -> meta
        }

        withProgress(message = 'beep boop, doing extraction things', value = 1, {
          test <- extract.metagene(
            as.character(meta[[6]]$genes),
            as.numeric(meta[[6]]$weights),
            beta2m(temp.processed$betas),
            as.numeric(meta[[7]])
          )
          round(test, digits = 3)
        })}),
        options = list(
          pageLength = 100,
          initComplete = I("function(settings, json) {alert('Done.');}"),
          processing=FALSE),
          #selection = list(target = 'row+column'),
          selection = "single")
        # test <- input$Mval_rows_selected
        


      
        output$time <- renderText({proc.time() - ptm})
        # cat(test)
        # test <- input$Mval_rows_all
        # test <- input$Mval_cell_clicked
        # s1 <- test$col
        # s2 <- test$value
        output$figure <-
          renderPlot(
            # test <- input$Mval_rows_all,

            # https://rstudio.github.io/DT/shiny.html

            generate_figure
                                    (c(#input$Mval_cells_selected,
                                       #"test$col" = test,
                                      #input$Mval_cell_clicked$value,
                                      "hi" = 1,
                                      "sup" = 0.4))
                                    # (c("patientA" = 0.1,
                                    #    test,
                                    #    "fljkfd" = 1))
        )
      from <- "Server"
      message <- "Processing complete!"
      messageData <- data.frame(from, message)
      output$messageMenu <- renderMenu({
        # Code to generate each of the messageItems here, in a list. This assumes
        # that messageData is a data frame with two columns, 'from' and 'message'.
        msgs <- apply(messageData, 1, function(row) {
          messageItem(from = row[["from"]], message = row[["message"]])
        })
        
        # This is equivalent to calling:
        #   dropdownMenu(type="messages", msgs[[1]], msgs[[2]], ...)
        dropdownMenu(type = "messages", badgeStatus = "success", .list = msgs)
      })
    })
  # output$figure <- renderPlot(generate_figure
  #                             #(input$input$tableId_row_last_clicked)
  #                             (c("patientA" = 0.1,
  #                                input$input$tableId_row_last_clicked,
  #                                "fljkfd" = 1))
  # )
  #Reset session and delete tempDIR
  observeEvent(input$bttn2, {
    #unlink(tempDIR, recursive = T)
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
