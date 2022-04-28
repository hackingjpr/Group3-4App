source("./source_functions.R")


ui <- shiny::fluidPage( 
  tags$h1("NAME OF APP TO FOLLOW"),
  
  setBackgroundColor(
    color = c("#F7FBFF", "#FFA500"),
    gradient = "radial",
    direction = c("top", "left")
  ),
  
  fluidRow(fileInput("idatFile", "Upload idat file:",multiple = TRUE, accept = ".idat"),
           
           column(3,actionBttn(
             inputId = "bttn1",
             label = "M value me!",
             color = "royal",
             style = "jelly"
           )),
           selectInput("metagenes", "Metagenes", c("ALL", "ATRT", "ECRT")),
           # column(3, tags$h2("Download"),
           #        radioButtons(inputId = "download", label = "Select file type", choices = c("png", "pdf")),
           #        downloadButton("down", "Download the results"),
                 # )),
           fluidRow( textOutput(outputId = "Mval")),
  tags$footer("Author: James Hacking,
                             Date Created: ?-?-2022,
                             Copyright (c) James Hacking, 2022,
              Email: james.hacking@ncl.ac.uk")))

server <- function(input, output){
  options(shiny.maxRequestSize=30*1024^2)
  observeEvent(input$bttn1,{  
    req(input$idatFile)
    

    
  output$Mval <- renderText({ 
    temp.base <- get_basenames(input$idatFile)
    # temp.processed <- process_idats(temp.base)
    # extract.metagene(as.character(atrt.meth.os.meta.n8.extract[[6]]$genes),
    # as.numeric(atrt.meth.os.meta.n8.extract[[6]]$weights),
    # beta2m(temp.processed$betas),
    # as.numeric(atrt.meth.os.meta.n8.extract[[7]]) )
  })})
}

shinyApp(server = server, ui = ui)



