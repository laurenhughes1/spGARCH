# -----------------------------
# Shiny app: Purely Spatial (No Time) GARCH/ARCH Snapshot
# - W built with base::toeplitz(), row-normalized
# - Pre-simulated epsilon^2 (and sigma^2) white noise ONCE at startup
# - Sliders: omega, alpha, beta (unused), gamma, lambda, rho (W decay), d = sqrt(N)
# - Outputs: ggplot heatmap and contour
# -----------------------------

library(shiny)
library(ggplot2)

# ---------- Helpers ----------

# Build Toeplitz W with exponential decay and row normalization
make_toeplitz <- function(N, rho = 0.5, normalize = TRUE) {
  first_col <- rho^(0:(N - 1))
  W <- toeplitz(first_col)
  if (normalize) {
    W <- W / rowSums(W)
  }
  W
}

# Compute sigma^2 from the purely spatial equation:
# (I - lambda W) sigma^2 = omega*1 + alpha*p + gamma*W p  where p = eps^2
compute_sigma2_spatial <- function(eps_vec, omega, alpha, gamma, W, lambda) {
  N <- length(eps_vec)
  p <- eps_vec^2
  rhs <- rep(omega, N) + alpha * p + gamma * as.vector(W %*% p)
  if (abs(lambda) < 1e-12) {
    return(rhs)
  } else {
    A <- diag(N) - lambda * W
    return(as.vector(solve(A, rhs)))
  }
}

# Make grid dataframe for ggplot (reshape a vector into d x d grid)
vec_to_grid_df <- function(v, d) {
  stopifnot(length(v) == d * d)
  df <- expand.grid(y = 1:d, x = 1:d)
  # fill by row so that increasing index fills across rows
  df$value <- as.vector(matrix(v, nrow = d, byrow = TRUE))
  df
}

# ---------- App Config ----------
# We allow d up to this; we will pre-simulate for max_N and take the first N = d^2 each time
D_MAX <- 40
MAX_N <- D_MAX * D_MAX
set.seed(123)

# Pre-simulate once (white noise). We store both epsilon and sigma2 noise;
# only epsilon is used in purely spatial model, sigma2_wn kept for completeness.
EPS_WN <- rnorm(MAX_N)
SIGMA2_WN <- abs(rnorm(MAX_N, mean = 1, sd = 0.2))^2

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("Purely Spatial GARCH/ARCH (No Time) — Toeplitz W"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("d", "Grid size d (so N = d\u00D7d):", min = 5, max = D_MAX, value = 20, step = 1),
      sliderInput("omega", "ω (scalar, replicated across units):", min = 0, max = 1, value = 0.1, step = 0.01),
      sliderInput("alpha", "α:", min = 0, max = 1, value = 0.2, step = 0.01),
      sliderInput("beta",  "β (unused here; temporal parameter):", min = 0, max = 0.99, value = 0.6, step = 0.01),
      sliderInput("gamma", "γ (weight on W p):", min = 0, max = 1, value = 0.1, step = 0.01),
      sliderInput("lambda","λ (spatial feedback in σ²):", min = 0, max = 0.95, value = 0.3, step = 0.01),
      sliderInput("rho",   "ρ (Toeplitz decay in W):", min = 0, max = 0.99, value = 0.5, step = 0.01),
      checkboxInput("normalizeW", "Row-normalize W", value = TRUE),
      helpText("Note: β is shown for completeness but not used (no temporal dynamics).")
    ),
    mainPanel(
      fluidRow(
        column(6, plotOutput("heatmapPlot", height = "450px")),
        column(6, plotOutput("contourPlot", height = "450px"))
      ),
      br(),
      verbatimTextOutput("stabilityInfo")
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  
  # Reactive: N and vectors sliced from pre-simulated pools (no re-simulation)
  N <- reactive({ input$d * input$d })
  
  eps_vec <- reactive({
    # Take first N entries from the global pre-simulated white noise
    EPS_WN[seq_len(N())]
  })
  
  # Build W (Toeplitz with row-normalization option)
  W <- reactive({
    make_toeplitz(N = N(), rho = input$rho, normalize = input$normalizeW)
  })
  
  # Compute sigma^2 snapshot
  sigma2_vec <- reactive({
    # Basic stability guard for (I - λW): ensure spectral radius(λ W) < 1
    eig <- eigen(W(), only.values = TRUE)$values
    rhoW <- max(Mod(eig))
    if (input$lambda * rhoW >= 1) {
      # Still compute (may error); we also warn in stability panel
      # To avoid crashing, nudge lambda slightly smaller than 1/rhoW
      lambda_safe <- min(input$lambda, 0.999 / rhoW)
      compute_sigma2_spatial(
        eps_vec(), input$omega, input$alpha, input$gamma, W(), lambda_safe
      )
    } else {
      compute_sigma2_spatial(
        eps_vec(), input$omega, input$alpha, input$gamma, W(), input$lambda
      )
    }
  })
  
  # Grid data for ggplot
  grid_df <- reactive({
    vec_to_grid_df(sigma2_vec(), d = input$d)
  })
  
  # Heatmap (ggplot)
  output$heatmapPlot <- renderPlot({
    df <- grid_df()
    ggplot(df, aes(x = x, y = y, fill = value)) +
      geom_tile() +
      coord_equal() +
      scale_fill_viridis_c(name = expression(sigma^2)) +
      labs(title = expression("Heatmap of " * sigma^2),
           x = "Spatial X", y = "Spatial Y") +
      theme_minimal(base_size = 13)
  })
  
  # Contour (ggplot)
  output$contourPlot <- renderPlot({
    df <- grid_df()
    ggplot(df, aes(x = x, y = y, z = value)) +
      geom_contour_filled(bins = 12, alpha = 0.9) +
      coord_equal() +
      labs(title = expression("Contour of " * sigma^2),
           x = "Spatial X", y = "Spatial Y", fill = expression(sigma^2)) +
      theme_minimal(base_size = 13)
  })
  
  # Stability info: show spectral radius check and reminders
  output$stabilityInfo <- renderPrint({
    eig <- eigen(W(), only.values = TRUE)$values
    rhoW <- max(Mod(eig))
    cat("N =", N(), "  (d =", input$d, ")\n")
    cat("spectral radius rho(W) ≈", round(rhoW, 4), "\n")
    cat("λ * rho(W) =", round(input$lambda * rhoW, 4),
        ifelse(input$lambda * rhoW < 1, "(OK: invertible I - λW)", "(WARNING: near/non-invertible)"), "\n")
    cat("Note: β is not used (no temporal dynamics).\n")
    invisible()
  })
}

shinyApp(ui, server)
