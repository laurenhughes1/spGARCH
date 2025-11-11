# ============================================================
# Shiny: SpGARCH(1,1) with Toeplitz W 
# - No re-simulation when parameters change
# - Pre-simulated ε (N_max x T_max) and initial σ^2 (N_max) at startup
# - Model:
#   σ_t^2 = ω·1 + α p_{t-1} + β σ_{t-1}^2 + γ W ε_{t-1}^2
# - Plots : heatmap + contour at selected time point
# - Constraint enforced: α + β + γ·ρ(W) < 1
# ============================================================

library(shiny)
library(ggplot2)

# ---------- Helpers ----------
# Geometric-decay Toeplitz W with diagonal = 0; optional row-normalization
make_toeplitz <- function(N, rho = 0.5, normalize = FALSE) {
  first_col <- rho^(0:(N - 1))
  W <- toeplitz(first_col)
  diag(W) <- 0  # enforce zero diagonal
  if (normalize) {
    rs <- rowSums(W)
    rs[rs == 0] <- 1  # avoid division-by-zero (e.g., when rho = 0)
    W <- W / rs
  }
  W
}

# Spectral radius ρ(W) (largest |eigenvalue|)
spectral_radius <- function(W) {
  max(abs(eigen(W, only.values = TRUE)$values))
}

# Reshape vector of length d*d into grid data frame for ggplot
vec_to_grid_df <- function(v, d) {
  stopifnot(length(v) == d * d)
  df <- expand.grid(y = 1:d, x = 1:d)
  df$value <- as.vector(matrix(v, nrow = d, byrow = TRUE))
  df
}

# ---------- App Config ----------
D_MAX  <- 40       # max d (so max N = D_MAX^2)
T_MAX  <- 1000     # max time to pre-simulate
set.seed(5)

# Pre-simulate once (global, reused) -- these are ε (not z)
EPS_POOL <- matrix(rnorm(D_MAX * D_MAX * T_MAX),
                   nrow = D_MAX * D_MAX, ncol = T_MAX)
SIGMA2_INIT_POOL <- abs(rnorm(D_MAX * D_MAX, mean = 1, sd = 0.25))^2

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("SpGARCH Simulator"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("d", "Grid size d (N = d×d):",
                  min = 5, max = D_MAX, value = 20, step = 1),
      numericInput("T", "T (number of time steps):",
                   value = 300, min = 5, max = T_MAX, step = 1),
      sliderInput("tstar", "Time point to visualize:",
                  min = 1, max = 300, value = 150, step = 1),
      hr(),
      sliderInput("omega", "ω (scalar):", min = 0, max = 1, value = 0.1, step = 0.01),
      sliderInput("alpha", "α:",          min = 0, max = 1, value = 0.2, step = 0.01),
      sliderInput("beta",  "β:",          min = 0, max = 0.99, value = 0.6, step = 0.01),
      sliderInput("gamma", "γ:",          min = 0, max = 1, value = 0.1, step = 0.01),
      hr(),
      sliderInput("rho",   "ρ (Toeplitz decay in W):",
                  min = 0, max = 0.99, value = 0.5, step = 0.01),
      checkboxInput("normalizeW", "Row-normalize W", value = FALSE),
      hr(),
      verbatimTextOutput("stabilityInfo")
    ),
    mainPanel(
      fluidRow(
        column(6, plotOutput("heatmapPlot", height = "450px")),
        column(6, plotOutput("contourPlot", height = "450px"))
      ),
      br()
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  
  # Keep the time slider in sync with T
  observeEvent(input$T, {
    newT <- max(5, min(input$T, T_MAX))
    updateSliderInput(session, "tstar",
                      min = 1, max = newT,
                      value = min(input$tstar, newT))
  })
  
  N <- reactive({ input$d * input$d })
  
  # Slice pre-simulated pools (no re-sim)
  eps_mat <- reactive({
    EPS_POOL[seq_len(N()), seq_len(input$T), drop = FALSE]
  })
  
  sigma2_init <- reactive({
    SIGMA2_INIT_POOL[seq_len(N())]
  })
  
  # Build W (geometric Toeplitz, diag=0)
  W <- reactive({
    make_toeplitz(N = N(), rho = input$rho, normalize = input$normalizeW)
  })
  
  # Stability constraint checker
  stability <- reactive({
    Wm <- W()
    rhoW <- spectral_radius(Wm)
    lhs  <- input$alpha + input$beta + input$gamma * rhoW
    list(rhoW = rhoW, lhs = lhs, ok = (lhs < 1 - 1e-10))
  })
  
  # Run temporal recursion deterministically with fixed eps and initial sigma2
  sigma2_cube <- reactive({
    # Enforce α + β + γ·ρ(W) < 1 (block rendering if violated)
    stab <- stability()
    validate(need(stab$ok,
                  sprintf("Stability constraint violated: α + β + γ·ρ(W) = %.4f (must be < 1). ρ(W)=%.4f",
                          stab$lhs, stab$rhoW)))
    
    Nn <- N(); Tt <- input$T
    Wm <- W()
    omega <- input$omega
    alpha <- input$alpha
    beta  <- input$beta
    gamma <- input$gamma
    
    sigma2 <- matrix(NA_real_, nrow = Nn, ncol = Tt)
    eps    <- eps_mat()
    
    # Initialize at t=1 with pre-simulated sigma2
    sigma2[, 1] <- sigma2_init()
    
    for (t in 2:Tt) {
      p_prev <- eps[, t - 1]^2
      sigma2[, t] <- rep(omega, Nn) +
        alpha * p_prev +
        beta  * sigma2[, t - 1] +
        gamma * as.vector(Wm %*% p_prev)
    }
    sigma2
  })
  
  # Data for plots at selected time point
  grid_df <- reactive({
    d <- input$d
    v <- sigma2_cube()[, input$tstar]
    vec_to_grid_df(v, d)
  })
  
  # Heatmap
  output$heatmapPlot <- renderPlot({
    df <- grid_df()
    ggplot(df, aes(x = x, y = y, fill = value)) +
      geom_tile() +
      coord_equal() +
      scale_fill_viridis_c(name = expression(sigma^2)) +
      labs(title = bquote("Heatmap of " * sigma^2 * " at t=" * .(input$tstar)),
           x = "Spatial X", y = "Spatial Y") +
      theme_minimal(base_size = 13)
  })
  
  # Contour
  output$contourPlot <- renderPlot({
    df <- grid_df()
    ggplot(df, aes(x = x, y = y, z = value)) +
      geom_contour_filled(bins = 12, alpha = 0.9) +
      coord_equal() +
      labs(title = bquote("Contour of " * sigma^2 * " at t=" * .(input$tstar)),
           x = "Spatial X", y = "Spatial Y", fill = expression(sigma^2)) +
      theme_minimal(base_size = 11)
  })
  
  # Diagnostics
  output$stabilityInfo <- renderPrint({
    stab <- stability()
    cat("N =", N(), "(d =", input$d, "),  T =", input$T, ",  t* =", input$tstar, "\n")
    cat("α + β =", round(input$alpha + input$beta, 4), "\n")
    cat(sprintf("ρ(W) = %.6f\n", stab$rhoW))
    cat(sprintf("α + β + γ·ρ(W) = %.6f  -> %s\n",
                stab$lhs, if (stab$ok) "OK (< 1)" else "VIOLATES (< 1 required)"))
    if (input$normalizeW) {
      cat("W is row-normalized (rows sum to 1).\n")
    } else {
      cat("W is NOT row-normalized.\n")
    }
    invisible()
  })
}

shinyApp(ui, server)
