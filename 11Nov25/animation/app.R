# ============================================================
# SpGARCH(1,1) Simulator — NO-JSON + Animation
# - Pre-simulated z ~ N(0,1)
# - eps = sigma * z
# - Toeplitz W (diag = 0) with rho slider
# - Stability: alpha + beta + gamma * rho(W) < 1
# - Animation: Play/Pause/Reset, FPS, Loop
# ============================================================

library(shiny)
library(ggplot2)

# ---------- Helpers ----------
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

# ---------- App Config ----------
D_MAX  <- 40
T_MAX  <- 1000
set.seed(5)

# Pre-simulated z (standard normals) and initial sigma^2
Z_POOL <- matrix(rnorm(D_MAX * D_MAX * T_MAX),
                 nrow = D_MAX * D_MAX, ncol = T_MAX)
SIGMA2_INIT_POOL <- abs(rnorm(D_MAX * D_MAX, mean = 1, sd = 0.25))^2

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("SpGARCH(1,1) Simulator — NO-JSON + Animation"),
  sidebarLayout(
    sidebarPanel(
      # Model size / time
      sliderInput("d", "Grid size d (N = d×d):",
                  min = 5, max = D_MAX, value = 20, step = 1),
      numericInput("T", "Number of time steps (T):",
                   value = 300, min = 5, max = T_MAX, step = 1),
      sliderInput("tstar", "Time point (t*):",
                  min = 1, max = 300, value = 150, step = 1),
      hr(),
      
      # Parameters
      sliderInput("omega", "ω (scalar):", min = 0, max = 1, value = 0.1, step = 0.01),
      sliderInput("alpha", "α:", min = 0, max = 1, value = 0.2, step = 0.01),
      sliderInput("beta",  "β:", min = 0, max = 0.99, value = 0.6, step = 0.01),
      sliderInput("gamma", "γ:", min = 0, max = 1, value = 0.1, step = 0.01),
      hr(),
      
      # W options
      sliderInput("rhoW", "ρ (Toeplitz decay in W):",
                  min = 0, max = 0.99, value = 0.5, step = 0.01),
      checkboxInput("normalizeW", "Row-normalize W", value = FALSE),
      hr(),
      
      # Plot & animation controls
      radioButtons("plot_var", "Quantity to show:",
                   choices = c("σ² (variance)" = "sigma2",
                               "ε (residuals)" = "eps"),
                   selected = "sigma2"),
      
      fluidRow(
        column(4, actionButton("play_pause", "▶ Play")),
        column(4, actionButton("reset_t", "⟲ Reset")),
        column(4, checkboxInput("loop", "Loop", value = TRUE))
      ),
      sliderInput("fps", "FPS (frames/sec):", min = 1, max = 30, value = 8, step = 1),
      
      hr(),
      verbatimTextOutput("stabilityInfo")
    ),
    
    mainPanel(
      fluidRow(
        column(6, plotOutput("heatmapPlot", height = "460px")),
        column(6, plotOutput("contourPlot", height = "460px"))
      )
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  
  # Keep t* bounds synced with T
  observeEvent(input$T, {
    updateSliderInput(session, "tstar",
                      min = 1, max = input$T,
                      value = min(input$tstar, input$T))
  })
  
  # Animation state
  playing <- reactiveVal(FALSE)
  
  observeEvent(input$play_pause, {
    playing(!isTRUE(playing()))
    updateActionButton(session, "play_pause", label = if (isTRUE(playing())) "⏸ Pause" else "▶ Play")
  })
  
  observeEvent(input$reset_t, {
    updateSliderInput(session, "tstar", value = 1)
  })
  
  # Frame clock
  observe({
    req(isTRUE(playing()))
    # Convert FPS to milliseconds per frame, clamp to reasonable min
    ms <- max(10, round(1000 / max(1, as.numeric(input$fps))))
    invalidateLater(ms, session)
    # advance t*
    t <- isolate(input$tstar)
    Tt <- isolate(input$T)
    if (t < Tt) {
      updateSliderInput(session, "tstar", value = t + 1)
    } else if (isTRUE(isolate(input$loop))) {
      updateSliderInput(session, "tstar", value = 1)
    } else {
      # stop at end if not looping
      playing(FALSE)
      updateActionButton(session, "play_pause", label = "▶ Play")
    }
  })
  
  # Core reactives
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
    validate(need(stab$ok,
                  sprintf("Stability violated: α+β+γ·ρ(W)=%.4f (must be <1). ρ(W)=%.4f",
                          stab$lhs, stab$rhoW)))
    
    Nn <- N(); Tt <- input$T
    Wm <- W()
    a <- input$alpha; b <- input$beta; g <- input$gamma; w <- input$omega
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
  
  # Grid data at current t*
  grid_df <- reactive({
    d <- input$d
    X <- if (input$plot_var == "sigma2") sim()$sigma2[, input$tstar] else sim()$eps[, input$tstar]
    vec_to_grid_df(X, d)
  })
  
  # Plots
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
  
  # Diagnostics
  output$stabilityInfo <- renderPrint({
    stab <- stability()
    cat(sprintf("ρ(W) = %.6f\n", stab$rhoW))
    cat(sprintf("α + β + γ·ρ(W) = %.6f  -> %s\n",
                stab$lhs, if (stab$ok) "OK (<1)" else "VIOLATES (<1 required)"))
    cat(sprintf("Playing: %s | t* = %d / %d | FPS = %d | Loop = %s\n",
                if (isTRUE(playing())) "Yes" else "No",
                input$tstar, input$T, as.integer(input$fps),
                if (isTRUE(input$loop)) "Yes" else "No"))
    invisible()
  })
}

shinyApp(ui, server)
