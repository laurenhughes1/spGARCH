# app.R
# Minimal Shiny JSON app using jsonlite::fromJSON(PATH) + diagnostics
# install.packages(c("shiny","jsonlite","DT")) if needed

library(shiny)
library(DT)

ui <- fluidPage(
  titlePanel("JSON Loader (Minimal)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("json", "Upload a JSON file",
                accept = c(".json", "application/json")),
      actionButton("loadSample", "Load sample JSON"),
      tags$hr(),
      downloadButton("downloadCsv", "Download parsed CSV")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Table", DTOutput("table")),
        tabPanel("Summary", verbatimTextOutput("summary")),
        tabPanel("Diagnostics", verbatimTextOutput("diag")),
        tabPanel("Raw (first 300 chars)", verbatimTextOutput("raw"))
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(use_sample = FALSE, sample_path = NULL)
  
  observeEvent(input$loadSample, {
    tf <- tempfile(fileext = ".json")
    writeLines('[{"id":1,"name":"Alice","age":31},
                 {"id":2,"name":"Bob","age":27},
                 {"id":3,"name":"Cara","age":35}]',
               tf, useBytes = TRUE)
    rv$sample_path <- tf
    rv$use_sample  <- TRUE
  })
  
  json_path <- reactive({
    if (rv$use_sample && !is.null(rv$sample_path)) {
      rv$sample_path
    } else {
      req(input$json)
      validate(need(nrow(input$json) == 1, "Please upload exactly one file."))
      input$json$datapath
    }
  })
  
  # --- Parse via file path (jsonlite accepts a path directly)
  parsed_df <- reactive({
    path <- json_path()
    obj <- tryCatch(
      jsonlite::fromJSON(path, simplifyDataFrame = TRUE),
      error = function(e) {
        validate(need(FALSE, paste("JSON parse error:", e$message)))
        NULL
      }
    )
    if (is.data.frame(obj)) {
      obj
    } else if (is.list(obj)) {
      # try to coerce common list-of-equal-length fields
      tryCatch(as.data.frame(obj), error = function(e) NULL)
    } else {
      NULL
    } %>%
      { validate(need(!is.null(.) && nrow(.) >= 1,
                      "Could not turn JSON into a table (expect an array of objects).")); . }
  })
  
  output$table <- renderDT({
    datatable(parsed_df(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$summary <- renderPrint({
    summary(parsed_df())
  })
  
  output$downloadCsv <- downloadHandler(
    filename = function() paste0("parsed-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv"),
    content  = function(file) write.csv(parsed_df(), file, row.names = FALSE)
  )
  
  output$raw <- renderText({
    path <- json_path()
    txt <- tryCatch(
      readChar(path, file.info(path)$size, useBytes = TRUE),
      error = function(e) paste("<could not read file:", e$message, ">")
    )
    substr(txt, 1, 300)
  })
  
  output$diag <- renderPrint({
    cat("Using path:\n", json_path(), "\n\n")
    cat("find('fromJSON'):\n"); print(find("fromJSON"))
    cat("\nSession packages:\n"); print(search())
  })
}

shinyApp(ui, server)
