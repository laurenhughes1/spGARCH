# app.R
if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("spdep", quietly = TRUE)) install.packages("spdep")
# install.packages("spGARCH")  # if needed

library(shiny)
library(ggplot2)
library(reshape2)
library(spdep)

make_W <- function(n_side, type = c("rook","queen")) {
  type <- match.arg(type)
  nb <- spdep::cell2nb(n_side, n_side, type = type)
  spdep::nb2mat(nb, style = "W")
}

ui <- fluidPage(
  titlePanel("Spatial GARCH simulation (spGARCH::sim.spGARCH)"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("n_side",  "Lattice size (n × n)", min = 10, max = 80, value = 40, step = 2),
      sliderInput("t_dummy", "Time steps (display only)", min = 1, max = 1, value = 1),  # sim.spGARCH returns 1 field
      tags$hr(),
      sliderInput("rho",    "ρ (rho) — ARCH spatial weight", min = -0.99, max = 0.99, value = 0.3, step = 0.01),
      sliderInput("lambda", "λ (lambda) — GARCH spatial weight", min = 0, max = 0.99, value = 0.5, step = 0.01),
      sliderInput("alpha",  "α (alpha) — unconditional variance level", min = 0, max = 3, value = 1, step = 0.01),
      checkboxInput("same_W2", "Use same W for W2 (GARCH) as W1 (ARCH)", TRUE),
      selectInput("W_type", "Adjacency", c("rook","queen"), selected = "rook"),
      checkboxInput("center_scale", "Center & scale for plotting", TRUE),
      tags$hr(),
      actionButton("reseed", "Reseed white noise"),
      helpText(
        "sim.spGARCH does not accept custom innovations. ",
        "We instead fix the RNG seed in control=list(seed=...). ",
        "Changing ρ, λ, α reuses the same seed so ε is identical; ",
        "click Reseed to draw a new ε."
      ),
      tags$p("Install if needed:", code("install.packages('spGARCH')"))
    ),
    mainPanel(
      fluidRow(
        column(6, h4("Heatmap"),  plotOutput("heatmap", height = "480px")),
        column(6, h4("Contour"),  plotOutput("contour", height = "480px"))
      )
    )
  )
)

server <- function(input, output, session) {
  # Rebuild W on size/type changes
  W_r <- reactive({
    make_W(input$n_side, input$W_type)
  })
  
  # Keep one seed that defines the white-noise realization
  seed_r <- reactiveVal(sample.int(1e6, 1))
  observeEvent(input$reseed, seed_r(sample.int(1e6, 1)))
  
  # Simulation that reuses the same ε by holding the seed fixed
  sim_grid <- reactive({
    req(W_r())
    validate(need(requireNamespace("spGARCH", quietly = TRUE),
                  "Please install the 'spGARCH' package."))
    
    W1 <- W_r()
    W2 <- if (isTRUE(input$same_W2)) W1 else make_W(input$n_side, if (input$W_type=="rook") "queen" else "rook")
    n  <- nrow(W1)
    
    y <- spGARCH::sim.spGARCH(
      n       = n,
      rho     = input$rho,
      lambda  = input$lambda,
      alpha   = input$alpha,
      W1      = W1,
      W2      = W2,
      type    = "spGARCH",
      control = list(seed = seed_r())  # <- keep ε fixed across parameter changes
    )
    
    # reshape vector -> n_side x n_side grid
    G <- matrix(as.numeric(y)[seq_len(n)], nrow = input$n_side, ncol = input$n_side, byrow = TRUE)
    
    if (isTRUE(input$center_scale)) {
      m <- mean(G); s <- sd(as.vector(G))
      if (is.finite(s) && s > 0) G <- (G - m)/s else G <- G - m
    }
    G
  })
  
  output$heatmap <- renderPlot({
    G <- sim_grid()
    df <- melt(G); names(df) <- c("y","x","z")
    ggplot(df, aes(x=x, y=y, fill=z)) +
      geom_raster(interpolate = TRUE) +
      scale_fill_viridis_c(option = "C") +
      coord_fixed() +
      labs(x=NULL, y=NULL,
           title = sprintf("Heatmap (ρ=%.2f, λ=%.2f, α=%.2f)", input$rho, input$lambda, input$alpha),
           fill = "Value") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold"))
  })
  
  output$contour <- renderPlot({
    G <- sim_grid()
    df <- melt(G); names(df) <- c("y","x","z")
    ggplot(df, aes(x=x, y=y, z=z)) +
      geom_contour(bins = 12) +
      coord_fixed() +
      labs(x=NULL, y=NULL,
           title = sprintf("Contour (ρ=%.2f, λ=%.2f, α=%.2f)", input$rho, input$lambda, input$alpha)) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold"))
  })
}

shinyApp(ui, server)
