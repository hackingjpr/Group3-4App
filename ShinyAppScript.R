source("./source_functions.R")


ui <- shiny::fluidPage(
  tags$h1("NAME OF APP TO FOLLOW"),
  
  setBackgroundColor(
    color = c("#F7FBFF", "#1E90FF"),
    gradient = "radial",
    direction = c("top", "left")
  ),
  
  fluidRow(
    fileInput(
      "idatFile",
      "Upload idat file:",
      multiple = TRUE,
      accept = ".idat"
    ),
    
    column(
      3,
      actionBttn(
        inputId = "bttn1",
        label = "M value me!",
        color = "primary",
        style = "fill"
      )
    ),
    column(
      3,
      actionBttn(
        inputId = "bttn2",
        label = "RESULTS?!",
        color = "success",
        style = "stretch"
      )
    ),
    selectInput("metagenes", "Metagenes", c("ALL", "ATRT", "ECRT")),
    # column(3, tags$h2("Download"),
    #        radioButtons(inputId = "download", label = "Select file type", choices = c("png", "pdf")),
    #        downloadButton("down", "Download the results"),
    # )),
    verbatimTextOutput("Mval", placeholder = TRUE),
    tags$footer(
      "Author: James Hacking,
                             Date Created: 28-04-2022,
                             Copyright (c) James Hacking, 2022,
              Email: james.hacking@ncl.ac.uk"
      
    )
  )
)

server <- function(input, output) {
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  observeEvent(input$bttn1, {
    req(input$idatFile)
    
    output$Mval <- renderPrint ({
      withProgress(message = 'beep boop, doing basey things', value = 1, {
    X <- input$idatFile
    #temp.base <- get_basenames(X$datapath)})
    temp.base <- get_basenames("/home/njh264/Idats/Epic/")})
    withProgress(message = "beep boop, doing idat things, trust me it will finish I just haven't worked out how to link the loading bar to progress",value = 0,{
      # for (i in 1:150){
      # incProgress(amount = 1/150)
      #   Sys.sleep(1)}    
    temp.processed <- process_idats(temp.base)})
    # atrt.meth.os.meta.n8.extract <- "ATRT"
    # ecrt.meth.os.meta.n20.extract <- "ECRT"
    # all.meth.os.meta.n54.extract <- "ALL"
      withProgress(message = 'beep boop, doing extraction things', value = 1, {
    extract.metagene(
      as.character(input$metagenes[[6]]$genes),
      as.numeric(input$metagenes[[6]]$weights),
      beta2m(temp.processed$betas),
      as.numeric(input$metagenes[[7]])
    )
      })
  })
  })
  # observeEvent(input$bttn2, {
  #   output$Mval <- renderPrint({
  #     Mvals
  #   })
  # })
}

shinyApp(server = server, ui = ui)
