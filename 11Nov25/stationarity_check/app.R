# ============================================================
# SpGARCH(1,1) Simulator
# - No re-simulation when parameters change
# - Pre-simulated ε (N_max x T_max) and initial σ^2 (N_max) at startup
# - Model:
#   σ_t^2 = ω·1 + α p_{t-1} + β σ_{t-1}^2 + γ W ε_{t-1}^2
# - Plots : heatmap + contour at selected time point
# - Constraint enforced: α + β + γ·ρ(W) < 1
# ============================================================

library(shiny)
library(ggplot2)

make_toeplitz <- function(N, rho = 0.5, normalize = FALSE) {
  first_col <- rho^(0:(N - 1))
  W <- toeplitz(first_col)
  diag(W) <- 0
  if (normalize) {
    rs <- rowSums(W); rs[rs == 0] <- 1
    W <- W / rs
  }
  W
}

vec_to_grid_df <- function(v, d) {
  stopifnot(length(v) == d * d)
  df <- expand.grid(y = 1:d, x = 1:d)
  df$value <- as.vector(matrix(v, nrow = d, byrow = TRUE))
  df
}

spectral_radius <- function(W) max(abs(eigen(W, only.values = TRUE)$values))

D_MAX <- 40; T_MAX <- 1000
set.seed(5)
Z_POOL <- matrix(rnorm(D_MAX * D_MAX * T_MAX),
                 nrow = D_MAX * D_MAX, ncol = T_MAX)
SIGMA2_INIT_POOL <- abs(rnorm(D_MAX * D_MAX, mean = 1, sd = 0.25))^2

ui <- fluidPage(
  titlePanel("SpGARCH(1,1) Simulator — NO-JSON v2"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("d", "Grid size d (N = d×d):", min = 5, max = D_MAX, value = 20, step = 1),
      numericInput("T", "Number of time steps (T):", value = 300, min = 5, max = T_MAX, step = 1),
      sliderInput("tstar", "Time point to visualize:", min = 1, max = 300, value = 150, step = 1),
      hr(),
      sliderInput("omega", "ω (scalar):", min = 0, max = 1, value = 0.1, step = 0.01),
      sliderInput("alpha", "α:", min = 0, max = 1, value = 0.2, step = 0.01),
      sliderInput("beta",  "β:", min = 0, max = 0.99, value = 0.6, step = 0.01),
      sliderInput("gamma", "γ:", min = 0, max = 1, value = 0.1, step = 0.01),
      hr(),
      sliderInput("rhoW", "ρ (Toeplitz decay in W):", min = 0, max = 0.99, value = 0.5, step = 0.01),
      checkboxInput("normalizeW", "Row-normalize W", value = FALSE),
      hr(),
      radioButtons("plot_var", "Quantity to plot:",
                   choices = c("σ² (variance)" = "sigma2", "ε (residuals)" = "eps"),
                   selected = "sigma2"),
      verbatimTextOutput("stabilityInfo")
    ),
    mainPanel(
      fluidRow(
        column(6, plotOutput("heatmapPlot", height = "450px")),
        column(6, plotOutput("contourPlot", height = "450px"))
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$T, {
    updateSliderInput(session, "tstar", min = 1, max = input$T, value = min(input$tstar, input$T))
  })
  
  N  <- reactive({ input$d * input$d })
  z  <- reactive({ Z_POOL[seq_len(N()), seq_len(input$T), drop = FALSE] })
  s0 <- reactive({ SIGMA2_INIT_POOL[seq_len(N())] })
  W  <- reactive({ make_toeplitz(N = N(), rho = input$rhoW, normalize = isTRUE(input$normalizeW)) })
  
  stability <- reactive({
    rhoW <- spectral_radius(W())
    lhs  <- input$alpha + input$beta + input$gamma * rhoW
    list(rhoW = rhoW, lhs = lhs, ok = (lhs < 1 - 1e-10))
  })
  
  sim <- reactive({
    stab <- stability()
    validate(need(stab$ok, sprintf("Stability violated: α+β+γ·ρ(W)=%.4f (must be <1). ρ(W)=%.4f",
                                   stab$lhs, stab$rhoW)))
    Nn <- N(); Tt <- input$T
    Wm <- W(); a <- input$alpha; b <- input$beta; g <- input$gamma; w <- input$omega
    sigma2 <- matrix(NA_real_, Nn, Tt)
    eps    <- matrix(NA_real_, Nn, Tt)
    zt     <- z()
    
    sigma2[,1] <- s0()
    eps[,1]    <- sqrt(sigma2[,1]) * zt[,1]
    
    for (t in 2:Tt) {
      p <- eps[, t-1]^2
      sigma2[, t] <- w + a * p + b * sigma2[, t-1] + g * as.vector(Wm %*% p)
      eps[, t]    <- sqrt(sigma2[, t]) * zt[, t]
    }
    list(sigma2 = sigma2, eps = eps)
  })
  
  grid_df <- reactive({
    d <- input$d
    x <- if (input$plot_var == "sigma2") sim()$sigma2[, input$tstar] else sim()$eps[, input$tstar]
    vec_to_grid_df(x, d)
  })
  
  output$heatmapPlot <- renderPlot({
    df <- grid_df()
    ttl <- if (input$plot_var == "sigma2")
      bquote("Heatmap of " * sigma^2 * " at t=" * .(input$tstar))
    else
      bquote("Heatmap of " * epsilon * " at t=" * .(input$tstar))
    ggplot(df, aes(x = x, y = y, fill = value)) +
      geom_tile() + coord_equal() +
      scale_fill_viridis_c(name = if (input$plot_var == "sigma2") expression(sigma^2) else expression(epsilon)) +
      labs(title = ttl, x = "Spatial X", y = "Spatial Y") +
      theme_minimal(base_size = 13)
  })
  
  output$contourPlot <- renderPlot({
    df <- grid_df()
    ttl <- if (input$plot_var == "sigma2")
      bquote("Contour of " * sigma^2 * " at t=" * .(input$tstar))
    else
      bquote("Contour of " * epsilon * " at t=" * .(input$tstar))
    ggplot(df, aes(x = x, y = y, z = value)) +
      geom_contour_filled(bins = 12, alpha = 0.9) +
      coord_equal() +
      labs(title = ttl, x = "Spatial X", y = "Spatial Y",
           fill = if (input$plot_var == "sigma2") expression(sigma^2) else expression(epsilon)) +
      theme_minimal(base_size = 11)
  })
  
  output$stabilityInfo <- renderPrint({
    stab <- stability()
    cat(sprintf("ρ(W) = %.6f\n", stab$rhoW))
    cat(sprintf("α + β + γ·ρ(W) = %.6f  -> %s\n",
                stab$lhs, if (stab$ok) "OK (<1)" else "VIOLATES (<1 required)"))
  })
}

shinyApp(ui, server)
