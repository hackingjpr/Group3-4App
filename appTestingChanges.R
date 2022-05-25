source("./source_functions.R")


ui <- shiny::fluidPage(
  
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
                                            title = "Risk Values",
                                            status = "info",
                                            solidHeader = TRUE,
                                            DTOutput('Mval')
                                          ))),
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
                                            h3("Metagene Set:"),
                                            textOutput("metagenechoice"),
                                            h3("Sample Selected:"),
                                            textOutput("sample"),
                                            h3("Patient's Risk Percentile:"),
                                            textOutput("percentages")
                                          ))
                                        )),
                               
                               #Download Tab
                               tabPanel("Download",
                                        textInput("filename", "Please insert desired filename", "Risk Values"),
                                        radioButtons(inputId = "download", label = "Select file type", choices = c("csv", "pdf")),
                                        downloadButton("down", "Download the results"),
                               )
                  )
                )),
                #Font Selection            
                tags$head(tags$style(HTML('* {font-family: "Courier New"};')))),
  
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
      att$set(10, text = "Loading") #Start at 10% 
      att$auto(ms = 1600, value = 0.01) # automatically increment
      
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
      
      test.res.mrt <- extract.metagene(
        as.character(ALL[[1]]$genes),
        as.numeric(ALL[[1]]$weights),
        beta2m(temp.processed$betas),
        as.numeric(ALL[[2]])
      )
      
      round(test.res.mrt, digits = 3) -> test.res.mrt
      
      
      test.res.atrt <- extract.metagene(
        as.character(ATRT[[1]]$genes),
        as.numeric(ATRT[[1]]$weights),
        beta2m(temp.processed$betas),
        as.numeric(ATRT[[2]])
      )
      
      round(test.res.atrt, digits = 3) -> test.res.atrt
      
      
      test.res.ecrt <- extract.metagene(
        as.character(ECRT[[1]]$genes),
        as.numeric(ECRT[[1]]$weights),
        beta2m(temp.processed$betas),
        as.numeric(ECRT[[2]])
      )
      
      round(test.res.ecrt, digits = 3) -> test.res.ecrt

      
      #print(test.res)
      # if (input$metagenes == "MRT (ATRT & ECRT)") {
      #   test.res.mrt -> test.res
      # }
      # if (input$metagenes == "ATRT") {
      #   test.res.atrt -> test.res
      # }
      # if (input$metagenes == "ECRT") {
      #   test.res.ecrt -> test.res
      # }
     # test.res <-  reactive( {
     #      if (input$metagenes == "MRT (ATRT & ECRT)") {
     #        test.res.mrt -> test.res
     #      }
     #      if (input$metagenes == "ATRT") {
     #        test.res.atrt -> test.res
     #      }
     #      if (input$metagenes == "ECRT") {
     #        test.res.ecrt -> test.res
     #      }
     #   })
      # message("test.res done")
      # print(test.res)
      
      test.reactive <- reactive(({
        test.res.mrt <- extract.metagene(
          as.character(ALL[[1]]$genes),
          as.numeric(ALL[[1]]$weights),
          beta2m(temp.processed$betas),
          as.numeric(ALL[[2]])
        )
        
        round(test.res.mrt, digits = 3) -> test.res.mrt
        
        
        test.res.atrt <- extract.metagene(
          as.character(ATRT[[1]]$genes),
          as.numeric(ATRT[[1]]$weights),
          beta2m(temp.processed$betas),
          as.numeric(ATRT[[2]])
        )
        
        round(test.res.atrt, digits = 3) -> test.res.atrt
        
        
        test.res.ecrt <- extract.metagene(
          as.character(ECRT[[1]]$genes),
          as.numeric(ECRT[[1]]$weights),
          beta2m(temp.processed$betas),
          as.numeric(ECRT[[2]])
        )
        
        round(test.res.ecrt, digits = 3) -> test.res.ecrt
        
        if (input$metagenes == "MRT (ATRT & ECRT)") {
          test.res.mrt -> test.res
        }
        if (input$metagenes == "ATRT") {
          test.res.atrt -> test.res
        }
        if (input$metagenes == "ECRT") {
          test.res.ecrt -> test.res
        }
        
        test.res}))
      
      output$Mval <- renderDT (
        ({
          test.res.mrt <- extract.metagene(
            as.character(ALL[[1]]$genes),
            as.numeric(ALL[[1]]$weights),
            beta2m(temp.processed$betas),
            as.numeric(ALL[[2]])
          )
          
          round(test.res.mrt, digits = 3) -> test.res.mrt
          
          
          test.res.atrt <- extract.metagene(
            as.character(ATRT[[1]]$genes),
            as.numeric(ATRT[[1]]$weights),
            beta2m(temp.processed$betas),
            as.numeric(ATRT[[2]])
          )
          
          round(test.res.atrt, digits = 3) -> test.res.atrt
          
          
          test.res.ecrt <- extract.metagene(
            as.character(ECRT[[1]]$genes),
            as.numeric(ECRT[[1]]$weights),
            beta2m(temp.processed$betas),
            as.numeric(ECRT[[2]])
          )
          
          round(test.res.ecrt, digits = 3) -> test.res.ecrt
          
          if (input$metagenes == "MRT (ATRT & ECRT)") {
            test.res.mrt -> test.res
          }
          if (input$metagenes == "ATRT") {
            test.res.atrt -> test.res
          }
          if (input$metagenes == "ECRT") {
            test.res.ecrt -> test.res
          }
          
          test.res}),
      options = list(
        pageLength = 10, 
        processing=FALSE),
      selection = "single"
      )

      att$done(text = "Complete")
      
      
      # output$time <- renderText({proc.time() - ptm})
      rowSelect <- reactive({input$Mval_rows_selected})
      message("row select done")
      
      output$figure <-
        renderPlot({
          
          
            test.res.mrt <- extract.metagene(
              as.character(ALL[[1]]$genes),
              as.numeric(ALL[[1]]$weights),
              beta2m(temp.processed$betas),
              as.numeric(ALL[[2]])
            )
            
            round(test.res.mrt, digits = 3) -> test.res.mrt 
            
            
            test.res.atrt <- extract.metagene(
              as.character(ATRT[[1]]$genes),
              as.numeric(ATRT[[1]]$weights),
              beta2m(temp.processed$betas),
              as.numeric(ATRT[[2]])
            )
            
            round(test.res.atrt, digits = 3) -> test.res.atrt
            
            
            test.res.ecrt <- extract.metagene(
              as.character(ECRT[[1]]$genes),
              as.numeric(ECRT[[1]]$weights),
              beta2m(temp.processed$betas),
              as.numeric(ECRT[[2]])
            )
            
            round(test.res.ecrt, digits = 3) -> test.res.ecrt
            
            if (input$metagenes == "MRT (ATRT & ECRT)") {
              test.res.mrt$Risk_Value -> figure.input
              names(figure.input) <- rownames(test.res.mrt)
            }
            if (input$metagenes == "ATRT") {
              test.res.atrt$Risk_Value -> figure.input
              names(figure.input) <- rownames(test.res.atrt)
            }
            if (input$metagenes == "ECRT") {
              test.res.ecrt$Risk_Value -> figure.input
              names(figure.input) <- rownames(test.res.ecrt)
            }
            
            # print(figure.input)
          
          # generate_figure_highlight_mrt(
          #   figure.input
          #   ,1)
          
          # generate_figure_highlight_mrt(
          #   figure.input
          #   ,input$Mval_row_last_clicked)
          # 
          figure.output.mrt <- (
            (if (input$metagenes == "MRT (ATRT & ECRT)") {
            # figureFile <- "./mrt54.dist.rds"
            generate_figure_highlight_mrt(
              figure.input
              ,input$Mval_row_last_clicked)
          }))
          figure.output.atrt <- (
            (if(input$metagenes == "ATRT") {
            # figureFile <- "./atrt8.dist.rds"
            generate_figure_highlight_atrt(
              figure.input
              ,input$Mval_row_last_clicked)
          }))
          figure.output.ecrt <-( 
            (if(input$metagenes == "ECRT") {
            # figureFile <- "./ecrt20.dist.rds"
            generate_figure_highlight_ecrt(
              figure.input
              ,input$Mval_row_last_clicked)
          }))
          # }
          
          if (input$metagenes == "MRT (ATRT & ECRT)") {
            figure.output.mrt -> figure.output
          }
          if (input$metagenes == "ATRT") {
            figure.output.atrt -> figure.output
          }
          if (input$metagenes == "ECRT") {
            figure.output.ecrt -> figure.output
          }
          
          
          figure.output
          # https://rstudio.github.io/DT/shiny.html
          
          
        })
      output$percentages <- renderText (
        {
        test.res.mrt <- extract.metagene(
          as.character(ALL[[1]]$genes),
          as.numeric(ALL[[1]]$weights),
          beta2m(temp.processed$betas),
          as.numeric(ALL[[2]])
        )
        
        round(test.res.mrt, digits = 3) -> test.res.mrt 
        
        
        test.res.atrt <- extract.metagene(
          as.character(ATRT[[1]]$genes),
          as.numeric(ATRT[[1]]$weights),
          beta2m(temp.processed$betas),
          as.numeric(ATRT[[2]])
        )
        
        round(test.res.atrt, digits = 3) -> test.res.atrt
        
        
        test.res.ecrt <- extract.metagene(
          as.character(ECRT[[1]]$genes),
          as.numeric(ECRT[[1]]$weights),
          beta2m(temp.processed$betas),
          as.numeric(ECRT[[2]])
        )
        
        round(test.res.ecrt, digits = 3) -> test.res.ecrt
        
        if (input$metagenes == "MRT (ATRT & ECRT)") {
          test.res.mrt$Risk_Value -> figure.input
          names(figure.input) <- rownames(test.res.mrt)
        }
        if (input$metagenes == "ATRT") {
          test.res.atrt$Risk_Value -> figure.input
          names(figure.input) <- rownames(test.res.atrt)
        }
        if (input$metagenes == "ECRT") {
          test.res.ecrt$Risk_Value -> figure.input
          names(figure.input) <- rownames(test.res.ecrt)
        }
        
        # print(figure.input)
        
        
        figure.output.mrt <- (
          (if (input$metagenes == "MRT (ATRT & ECRT)") {
            # figureFile <- "./mrt54.dist.rds"
            generate_figure_percentage_mrt(
              figure.input
              ,input$Mval_row_last_clicked)
          }))
        figure.output.atrt <- (
          (if(input$metagenes == "ATRT") {
            # figureFile <- "./atrt8.dist.rds"
            generate_figure_percentage_atrt(
              figure.input
              ,input$Mval_row_last_clicked)
          }))
        figure.output.ecrt <-( 
          (if(input$metagenes == "ECRT") {
            # figureFile <- "./ecrt20.dist.rds"
            generate_figure_percentage_ecrt(
              figure.input
              ,input$Mval_row_last_clicked)
          }))
        # }
        
        if (input$metagenes == "MRT (ATRT & ECRT)") {
          figure.output.mrt -> figure.output.percentage
        }
        if (input$metagenes == "ATRT") {
          figure.output.atrt -> figure.output.percentage
        }
        if (input$metagenes == "ECRT") {
          figure.output.ecrt -> figure.output.percentage
        }
        
        
        figure.output.percentage
        
        # input <- "./mrt54.dist.rds",
        # generate_figure_percentage_mrt(
        #   figure.input,
        #   input$Mval_row_last_clicked)
        # ,input$Mval_row_last_clicked)
        # "hello",
        
         # if (input$metagenes == "MRT (ATRT & ECRT)") {
         #   # figureFile <- "./mrt54.dist.rds"
         #   generate_figure_percentage_mrt(
         #     figure.input
         #     ,input$Mval_row_last_clicked)
         # }
         # if (input$metagenes == "ATRT") {
         #   # figureFile <- "./atrt8.dist.rds"
         #   generate_figure_percentage_atrt(
         #     figure.input
         #     ,input$Mval_row_last_clicked)
         # }
         # if (input$metagenes == "ECRT") {
         #   # figureFile <- "./ecrt20.dist.rds"
         #   generate_figure_percentage_ecrt(
         #     figure.input
         #     ,input$Mval_row_last_clicked)
         # }
      })

      output$metagenechoice <- renderText({input$metagenes})
      
      showTab(inputId = "tabs", target = "Results Table", select = TRUE)
      showTab(inputId = "tabs", target = "Download")
      
      #For displaying currently selected sample
      
      output$sample <- renderText(input$Mval_cell_clicked$value) 
      
      output$down <- downloadHandler(
        filename = function() {
          paste(input$filename,Sys.time(), input$download, sep=".")
        },
        content ={ 
          function(file) {
            if (input$download == "csv")
            write.csv(
              (
                ({
                  test.res.mrt <- extract.metagene(
                    as.character(ALL[[1]]$genes),
                    as.numeric(ALL[[1]]$weights),
                    beta2m(temp.processed$betas),
                    as.numeric(ALL[[2]])
                  )
                  
                  round(test.res.mrt, digits = 3) -> test.res.mrt
                  
                  
                  test.res.atrt <- extract.metagene(
                    as.character(ATRT[[1]]$genes),
                    as.numeric(ATRT[[1]]$weights),
                    beta2m(temp.processed$betas),
                    as.numeric(ATRT[[2]])
                  )
                  
                  round(test.res.atrt, digits = 3) -> test.res.atrt
                  
                  
                  test.res.ecrt <- extract.metagene(
                    as.character(ECRT[[1]]$genes),
                    as.numeric(ECRT[[1]]$weights),
                    beta2m(temp.processed$betas),
                    as.numeric(ECRT[[2]])
                  )
                  
                  round(test.res.ecrt, digits = 3) -> test.res.ecrt
                  
                  if (input$metagenes == "MRT (ATRT & ECRT)") {
                    test.res.mrt -> test.res
                  }
                  if (input$metagenes == "ATRT") {
                    test.res.atrt -> test.res
                  }
                  if (input$metagenes == "ECRT") {
                    test.res.ecrt -> test.res
                  }
                  
                  test.res}) 
                
              ), file, row.names = TRUE)
            
            else
              pdf(file,
                  width = 14
                  )
            grid.table(
              
              ({
                test.res.mrt <- extract.metagene(
                  as.character(ALL[[1]]$genes),
                  as.numeric(ALL[[1]]$weights),
                  beta2m(temp.processed$betas),
                  as.numeric(ALL[[2]])
                )
                
                round(test.res.mrt, digits = 3) -> test.res.mrt
                
                
                test.res.atrt <- extract.metagene(
                  as.character(ATRT[[1]]$genes),
                  as.numeric(ATRT[[1]]$weights),
                  beta2m(temp.processed$betas),
                  as.numeric(ATRT[[2]])
                )
                
                round(test.res.atrt, digits = 3) -> test.res.atrt
                
                
                test.res.ecrt <- extract.metagene(
                  as.character(ECRT[[1]]$genes),
                  as.numeric(ECRT[[1]]$weights),
                  beta2m(temp.processed$betas),
                  as.numeric(ECRT[[2]])
                )
                
                round(test.res.ecrt, digits = 3) -> test.res.ecrt
                
                if (input$metagenes == "MRT (ATRT & ECRT)") {
                  test.res.mrt -> test.res
                }
                if (input$metagenes == "ATRT") {
                  test.res.atrt -> test.res
                }
                if (input$metagenes == "ECRT") {
                  test.res.ecrt -> test.res
                }
                
                test.res})
            )
            plot({
              
              
                
                test.res.mrt <- extract.metagene(
                  as.character(ALL[[1]]$genes),
                  as.numeric(ALL[[1]]$weights),
                  beta2m(temp.processed$betas),
                  as.numeric(ALL[[2]])
                )
                
                round(test.res.mrt, digits = 3) -> test.res.mrt 
                
                
                test.res.atrt <- extract.metagene(
                  as.character(ATRT[[1]]$genes),
                  as.numeric(ATRT[[1]]$weights),
                  beta2m(temp.processed$betas),
                  as.numeric(ATRT[[2]])
                )
                
                round(test.res.atrt, digits = 3) -> test.res.atrt
                
                
                test.res.ecrt <- extract.metagene(
                  as.character(ECRT[[1]]$genes),
                  as.numeric(ECRT[[1]]$weights),
                  beta2m(temp.processed$betas),
                  as.numeric(ECRT[[2]])
                )
                
                round(test.res.ecrt, digits = 3) -> test.res.ecrt
                
                if (input$metagenes == "MRT (ATRT & ECRT)") {
                  test.res.mrt$Risk_Value -> figure.input
                  names(figure.input) <- rownames(test.res.mrt)
                }
                if (input$metagenes == "ATRT") {
                  test.res.atrt$Risk_Value -> figure.input
                  names(figure.input) <- rownames(test.res.atrt)
                }
                if (input$metagenes == "ECRT") {
                  test.res.ecrt$Risk_Value -> figure.input
                  names(figure.input) <- rownames(test.res.ecrt)
                }
                
                # print(figure.input)
                
                # generate_figure_highlight_mrt(
                #   figure.input
                #   ,1)
                
                # generate_figure_highlight_mrt(
                #   figure.input
                #   ,input$Mval_row_last_clicked)
                # 
                figure.output.mrt <- (
                  (if (input$metagenes == "MRT (ATRT & ECRT)") {
                    # figureFile <- "./mrt54.dist.rds"
                    generate_figure_highlight_mrt(
                      figure.input
                      ,input$Mval_row_last_clicked)
                  }))
                figure.output.atrt <- (
                  (if(input$metagenes == "ATRT") {
                    # figureFile <- "./atrt8.dist.rds"
                    generate_figure_highlight_atrt(
                      figure.input
                      ,input$Mval_row_last_clicked)
                  }))
                figure.output.ecrt <-( 
                  (if(input$metagenes == "ECRT") {
                    # figureFile <- "./ecrt20.dist.rds"
                    generate_figure_highlight_ecrt(
                      figure.input
                      ,input$Mval_row_last_clicked)
                  }))
                # }
                
                if (input$metagenes == "MRT (ATRT & ECRT)") {
                  figure.output.mrt -> figure.output
                }
                if (input$metagenes == "ATRT") {
                  figure.output.atrt -> figure.output
                }
                if (input$metagenes == "ECRT") {
                  figure.output.ecrt -> figure.output
                }
                
                
                figure.output
                
                
              })
              dev.off()
            
          }
          # if (input$download == "csv")
          #   write.csv(output$Mval)
          # else
          #   pdf(output$Mval,
          #       width = 14
          #   )
          # dev.off()
        })
      
      
    })
  #Reset session and delete tempDIR
  observeEvent(input$bttn2, {
    unlink(tempDIR, recursive = T)
    session$reload()
    return()
    print("session reload not working")
  })
  renderText(output$metagenes <- input$metagenes)
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