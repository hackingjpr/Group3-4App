if(!require(data.table)){
  install.packages("shiny")
  library(shiny)
}

if(!require(minfiData)){
  install.packages("minfiData")
  library(minfiData)
}

if(!require(sva)){
  install.packages("sva")
  library(sva)
}

if(!require(devtools)){
  install.packages("devtools")
  library(devtools)
}

if(!require(bumphunter)){
  install.packages("bumphunter")
  library(bumphunter)
}

if(!require(lumi)){
  install.packages("lumi")
  library(lumi)
}

if(!require(fluidPage)){
  install.packages("fluidPage")
  library(fluidPage)
}



ui <- shiny::fluidPage( 
  tags$h1("NAME OF APP TO FOLLOW"),
  fluidRow(fileInput("idatFile", "Upload idat file:", accept = ".idat")),
  tags$footer("Author: James Hacking,
                             Date Created: ?-?-2022,
                             Copyright (c) James Hacking, 2022,
              Email: james.hacking@ncl.ac.uk"))

server <- function(input, output)

}

shinyApp(server = server, ui = ui)



