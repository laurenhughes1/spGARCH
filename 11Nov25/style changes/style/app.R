# ============================================================
# SpGARCH(1,1) — Polished UI with Tabs, Presets, JSON & Animation
# - Default: Toeplitz W (diag=0), optional row-normalization
# - Optional: Upload W (JSON N×N), optional row-normalization
# - ω: scalar slider or upload vector (JSON length N=d^2)
# - Auto-adjust d from uploaded ω/W sizes (with validation badges)
# - z ~ N(0,1) pre-simulated; eps = sigma * z
# - Stability: alpha + beta + gamma * rho(W) < 1
# - Animation: Play/Pause/Reset, FPS, Loop
# - Nicer layout with navbar tabs & presets
# ============================================================

library(shiny)
library(ggplot2)
library(jsonlite)

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

safe_fromJSON <- function(path, simplifyVector = TRUE) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  txt <- try(readLines(path, warn = FALSE), silent = TRUE)
  if (inherits(txt, "try-error") || length(txt) == 0) return(NULL)
  out <- try(jsonlite::fromJSON(txt, simplifyVector = simplifyVector), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

# ---------- App Config ----------
D_MAX  <- 40      # pools preallocated to D_MAX^2
T_MAX  <- 1000
set.seed(5)

Z_POOL <- matrix(rnorm(D_MAX * D_MAX * T_MAX),
                 nrow = D_MAX * D_MAX, ncol = T_MAX)
SIGMA2_INIT_POOL <- abs(rnorm(D_MAX * D_MAX, mean = 1, sd = 0.25))^2

# ---------- UI ----------
ui <- navbarPage(
  title = "SpGARCH(1,1) — Spatial Volatility Lab",
  id = "nav",
  inverse = TRUE,  # dark header
  header = NULL,
  
  # --- SIMULATE (plots + most-used controls) ---
  tabPanel(
    "Simulate",
    fluidPage(
      br(),
      fluidRow(
        column(
          4,
          wellPanel(
            h4("Quick Presets"),
            selectInput("preset", label = NULL,
                        choices = c("Choose…" = "",
                                    "Small, fast (d=15, T=250, FPS=12, ρ=0.4)" = "small",
                                    "Medium (d=25, T=400, FPS=8, ρ=0.5)" = "medium",
                                    "Large, slow (d=35, T=600, FPS=6, ρ=0.6)" = "large")),
            actionButton("apply_preset", "Apply Preset", class = "btn-primary")
          ),
          wellPanel(
            h4("Layout & Time"),
            sliderInput("d", "Grid size d (N = d×d):",
                        min = 5, max = D_MAX, value = 20, step = 1),
            numericInput("T", "Number of time steps (T):",
                         value = 300, min = 5, max = T_MAX, step = 1),
            sliderInput("tstar", "Time point (t*):",
                        min = 1, max = 300, value = 150, step = 1)
          ),
          wellPanel(
            h4("What to Plot"),
            radioButtons("plot_var", label = NULL,
                         choices = c("σ² (variance)" = "sigma2", "ε (residuals)" = "eps"),
                         inline = TRUE, selected = "sigma2")
          ),
          wellPanel(
            h4("Animation"),
            fluidRow(
              column(4, actionButton("play_pause", "▶ Play", width = "100%")),
              column(4, actionButton("reset_t", "⟲ Reset", width = "100%")),
              column(4, checkboxInput("loop", "Loop", value = TRUE))
            ),
            sliderInput("fps", "FPS (frames/sec):", min = 1, max = 30, value = 8, step = 1)
          )
        ),
        column(
          8,
          fluidRow(
            column(12,
                   wellPanel(
                     h4("Plots"),
                     fluidRow(
                       column(6, plotOutput("heatmapPlot", height = "460px")),
                       column(6, plotOutput("contourPlot", height = "460px"))
                     )
                   )
            )
          )
        )
      )
    )
  ),
  
  # --- UPLOADS (ω & W JSON with badges) ---
  tabPanel(
    "Uploads",
    fluidPage(
      br(),
      fluidRow(
        column(
          6,
          wellPanel(
            h4("ω (Intercept)"),
            radioButtons("omega_src", "Source:", inline = TRUE,
                         choices = c("Scalar" = "scalar", "Upload JSON" = "upload_vec"),
                         selected = "scalar"),
            conditionalPanel(
              condition = "input.omega_src == 'scalar'",
              sliderInput("omega_scalar", "ω (scalar):", min = 0, max = 1, value = 0.1, step = 0.01)
            ),
            conditionalPanel(
              condition = "input.omega_src == 'upload_vec'",
              fileInput("omega_file", "Upload ω vector (JSON, length N=d^2)", accept = c(".json")),
              div(style = "margin-top: 6px;", uiOutput("omega_status"))
            )
          )
        ),
        column(
          6,
          wellPanel(
            h4("W (Spatial Weights)"),
            radioButtons("W_src", "Source:", inline = TRUE,
                         choices = c("Toeplitz (default)" = "toeplitz", "Upload JSON" = "upload_W"),
                         selected = "toeplitz"),
            conditionalPanel(
              condition = "input.W_src == 'toeplitz'",
              sliderInput("rhoW", "ρ (Toeplitz decay in W):",
                          min = 0, max = 0.99, value = 0.5, step = 0.01),
              checkboxInput("normalizeW", "Row-normalize Toeplitz W", value = FALSE),
              helpText("Toeplitz W uses diag=0.")
            ),
            conditionalPanel(
              condition = "input.W_src == 'upload_W'",
              fileInput("W_file", "Upload W matrix (JSON N×N numeric)", accept = c(".json")),
              checkboxInput("normalizeW_up", "Row-normalize uploaded W", value = FALSE),
              div(style = "margin-top: 6px;", uiOutput("W_status"))
            )
          )
        )
      ),
      fluidRow(
        column(
          12,
          wellPanel(
            h4("Auto grid size (d) from uploads"),
            helpText("When a valid ω (length N) and/or W (N×N) are uploaded, the app infers d=√N and updates the slider. If both disagree, no change is made and you’ll see a notification.")
          )
        )
      )
    )
  ),
  
  # --- MODEL (parameters & stability) ---
  tabPanel(
    "Model",
    fluidPage(
      br(),
      fluidRow(
        column(
          6,
          wellPanel(
            h4("GARCH Parameters"),
            sliderInput("alpha", "α", min = 0, max = 1, value = 0.2, step = 0.01),
            sliderInput("beta",  "β", min = 0, max = 0.99, value = 0.6, step = 0.01),
            sliderInput("gamma", "γ", min = 0, max = 1, value = 0.1, step = 0.01),
            helpText("Stability constraint enforced:  α + β + γ·ρ(W) < 1")
          )
        ),
        column(
          6,
          wellPanel(
            h4("Stability & Info"),
            verbatimTextOutput("stabilityInfo", placeholder = TRUE)
          )
        )
      )
    )
  ),
  
  # --- ANIMATION (focused controls) ---
  tabPanel(
    "Animation",
    fluidPage(
      br(),
      fluidRow(
        column(
          4,
          wellPanel(
            h4("Controls"),
            sliderInput("tstar", "Time (t*)", min = 1, max = 300, value = 150, step = 1),
            fluidRow(
              column(4, actionButton("play_pause", "▶ Play", width = "100%")),
              column(4, actionButton("reset_t", "⟲ Reset", width = "100%")),
              column(4, checkboxInput("loop", "Loop", value = TRUE))
            ),
            sliderInput("fps", "FPS", min = 1, max = 30, value = 8, step = 1)
          )
        ),
        column(
          8,
          wellPanel(
            h4("Live View"),
            plotOutput("heatmapPlot", height = "460px")
          )
        )
      )
    )
  ),
  
  # --- DIAGNOSTICS (plain text) ---
  tabPanel(
    "Diagnostics",
    fluidPage(
      br(),
      wellPanel(
        h4("Console"),
        verbatimTextOutput("stabilityInfo_diag", placeholder = TRUE)
      )
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  
  # ---- Presets ----
  observeEvent(input$apply_preset, {
    if (nzchar(input$preset)) {
      switch(input$preset,
             "small" = {
               updateSliderInput(session, "d", value = 15)
               updateNumericInput(session, "T", value = 250)
               updateSliderInput(session, "fps", value = 12)
               updateSliderInput(session, "rhoW", value = 0.4)
             },
             "medium" = {
               updateSliderInput(session, "d", value = 25)
               updateNumericInput(session, "T", value = 400)
               updateSliderInput(session, "fps", value = 8)
               updateSliderInput(session, "rhoW", value = 0.5)
             },
             "large" = {
               updateSliderInput(session, "d", value = 35)
               updateNumericInput(session, "T", value = 600)
               updateSliderInput(session, "fps", value = 6)
               updateSliderInput(session, "rhoW", value = 0.6)
             }
      )
      showNotification("Preset applied.", type = "message", duration = 2)
    }
  })
  
  # Sync t* with T
  observeEvent(input$T, {
    updateSliderInput(session, "tstar",
                      min = 1, max = input$T,
                      value = min(input$tstar, input$T))
  })
  
  # Animation controls
  playing <- reactiveVal(FALSE)
  observeEvent(input$play_pause, {
    playing(!isTRUE(playing()))
    updateActionButton(session, "play_pause",
                       label = if (isTRUE(playing())) "⏸ Pause" else "▶ Play")
  })
  observeEvent(input$reset_t, {
    updateSliderInput(session, "tstar", value = 1)
  })
  observe({
    req(isTRUE(playing()))
    ms <- max(10, round(1000 / max(1, as.numeric(input$fps))))
    invalidateLater(ms, session)
    t <- isolate(input$tstar); Tt <- isolate(input$T)
    if (t < Tt) updateSliderInput(session, "tstar", value = t + 1)
    else if (isTRUE(isolate(input$loop))) updateSliderInput(session, "tstar", value = 1)
    else {
      playing(FALSE)
      updateActionButton(session, "play_pause", label = "▶ Play")
    }
  })
  
  # Core sizes/data
  N  <- reactive({ input$d * input$d })
  z  <- reactive({ Z_POOL[seq_len(N()), seq_len(input$T), drop = FALSE] })
  s0 <- reactive({ SIGMA2_INIT_POOL[seq_len(N())] })
  
  # ---- ω: scalar or uploaded JSON (with badge) ----
  omega_parsed <- reactive({
    if (!identical(input$omega_src, "upload_vec")) return(NULL)
    if (is.null(input$omega_file) || !file.exists(input$omega_file$datapath)) return(NULL)
    vec <- safe_fromJSON(input$omega_file$datapath, simplifyVector = TRUE)
    suppressWarnings(as.numeric(vec))
  })
  
  output$omega_status <- renderUI({
    Nn <- N()
    if (!identical(input$omega_src, "upload_vec"))
      return(tags$span(class = "label label-default", "ω: Scalar mode"))
    if (is.null(input$omega_file))
      return(tags$span(class = "label label-default", "ω JSON: no file"))
    vec <- omega_parsed()
    if (is.null(vec) || any(!is.finite(vec)))
      return(tags$span(class = "label label-danger", "ω JSON: invalid"))
    if (length(vec) != Nn)
      return(tags$span(class = "label label-warning",
                       paste0("ω JSON length ", length(vec), " (current N=", Nn, ")")))
    tags$span(class = "label label-success", paste0("ω JSON matches N=", Nn))
  })
  
  omega_vec <- reactive({
    Nn <- N()
    omega_default <- rep(input$omega_scalar, Nn)
    if (!identical(input$omega_src, "upload_vec")) return(omega_default)
    vec <- omega_parsed()
    if (is.null(vec) || any(!is.finite(vec))) return(omega_default)
    if (length(vec) != Nn) return(omega_default)
    vec
  })
  
  # ---- W: Toeplitz or uploaded JSON (with badge) ----
  W_parsed <- reactive({
    if (!identical(input$W_src, "upload_W")) return(NULL)
    if (is.null(input$W_file) || !file.exists(input$W_file$datapath)) return(NULL)
    raw <- safe_fromJSON(input$W_file$datapath, simplifyVector = FALSE)
    if (is.null(raw)) return(NULL)
    Wm <- try({
      if (is.matrix(raw)) raw else
        do.call(rbind, lapply(raw, function(row) unlist(row, use.names = FALSE)))
    }, silent = TRUE)
    if (inherits(Wm, "try-error")) return(NULL)
    Wm <- as.matrix(Wm); storage.mode(Wm) <- "double"
    Wm
  })
  
  output$W_status <- renderUI({
    Nn <- N()
    if (!identical(input$W_src, "upload_W"))
      return(tags$span(class = "label label-default", "W: Toeplitz mode"))
    if (is.null(input$W_file))
      return(tags$span(class = "label label-default", "W JSON: no file"))
    Wm <- W_parsed()
    if (is.null(Wm) || any(!is.finite(Wm)))
      return(tags$span(class = "label label-danger", "W JSON: invalid"))
    tag <- if (nrow(Wm) == ncol(Wm)) "square" else "non-square"
    if (!(nrow(Wm) == Nn && ncol(Wm) == Nn))
      return(tags$span(class = "label label-warning",
                       paste0("W JSON ", tag, " ",
                              nrow(Wm), "×", ncol(Wm), " (current N=", Nn, ")")))
    tags$span(class = "label label-success", paste0("W JSON matches N×N (N=", Nn, ")"))
  })
  
  # ---- Auto d from uploads ----
  observeEvent(list(input$omega_file, input$W_file, input$omega_src, input$W_src), {
    d_from_omega <- NULL; d_from_W <- NULL
    if (identical(input$omega_src, "upload_vec")) {
      vec <- omega_parsed()
      if (!is.null(vec) && all(is.finite(vec))) {
        N_om <- length(vec); d_om <- sqrt(N_om)
        if (abs(d_om - round(d_om)) < 1e-9) d_from_omega <- as.integer(round(d_om))
      }
    }
    if (identical(input$W_src, "upload_W")) {
      Wm <- W_parsed()
      if (!is.null(Wm) && all(is.finite(Wm)) && nrow(Wm) == ncol(Wm)) {
        d_from_W <- as.integer(round(sqrt(nrow(Wm))))
        if (d_from_W * d_from_W != nrow(Wm)) d_from_W <- NULL
      }
    }
    d_target <- NULL
    if (!is.null(d_from_W) && !is.null(d_from_omega)) {
      if (d_from_W == d_from_omega) d_target <- d_from_W
      else {
        showNotification(
          paste0("ω and W sizes disagree (dω=", d_from_omega, ", dW=", d_from_W, "). Not changing d."),
          type = "error", duration = 6
        )
        return()
      }
    } else if (!is.null(d_from_W)) d_target <- d_from_W
    else if (!is.null(d_from_omega)) d_target <- d_from_omega
    else return()
    
    if (d_target > D_MAX) {
      showNotification(
        paste0("Uploaded data implies d=", d_target, " > D_MAX=", D_MAX,
               ". Increase D_MAX in code to use this size."),
        type = "error", duration = 6
      )
      return()
    }
    if (!identical(input$d, d_target)) {
      updateSliderInput(session, "d", min = 5, max = D_MAX, value = d_target)
      showNotification(paste0("Auto-set d to ", d_target, " to match uploaded data."), type = "message", duration = 4)
    }
  }, ignoreInit = TRUE)
  
  # Effective W
  W_effective <- reactive({
    Nn <- N()
    if (!identical(input$W_src, "upload_W")) {
      return(make_toeplitz(N = Nn, rho = input$rhoW, normalize = isTRUE(input$normalizeW)))
    }
    Wm <- W_parsed()
    if (is.null(Wm) || any(!is.finite(Wm)) || !(nrow(Wm) == Nn && ncol(Wm) == Nn)) {
      return(make_toeplitz(N = Nn, rho = input$rhoW, normalize = isTRUE(input$normalizeW)))
    }
    if (isTRUE(input$normalizeW_up)) {
      rs <- rowSums(Wm); rs[rs == 0] <- 1
      Wm <- Wm / rs
    }
    Wm
  })
  
  # Stability & simulation
  stability <- reactive({
    rhoW <- spectral_radius(W_effective())
    lhs  <- input$alpha + input$beta + input$gamma * rhoW
    list(rhoW = rhoW, lhs = lhs, ok = (lhs < 1 - 1e-10))
  })
  
  sim <- reactive({
    stab <- stability()
    validate(need(stab$ok,
                  sprintf("Stability violated: α+β+γ·ρ(W)=%.4f (must be <1). ρ(W)=%.4f",
                          stab$lhs, stab$rhoW)))
    
    Nn <- N(); Tt <- input$T
    Wm <- W_effective()
    a <- input$alpha; b <- input$beta; g <- input$gamma
    omegaV <- omega_vec()
    sigma2 <- matrix(NA_real_, Nn, Tt)
    eps    <- matrix(NA_real_, Nn, Tt)
    zt     <- z()
    
    sigma2[,1] <- s0()
    eps[,1]    <- sqrt(sigma2[,1]) * zt[,1]
    
    for (t in 2:Tt) {
      p <- eps[, t-1]^2
      sigma2[, t] <- omegaV + a * p + b * sigma2[, t-1] + g * as.vector(Wm %*% p)
      eps[, t]    <- sqrt(sigma2[, t]) * zt[, t]
    }
    list(sigma2 = sigma2, eps = eps)
  })
  
  # Grid & plots
  grid_df <- reactive({
    d <- input$d
    X <- if (input$plot_var == "sigma2") sim()$sigma2[, input$tstar] else sim()$eps[, input$tstar]
    vec_to_grid_df(X, d)
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
  
  # Diagnostics outputs in two places
  output$stabilityInfo <- renderPrint({
    stab <- stability()
    Nn <- N()
    cat(sprintf("N = %d  (d = %d),  T = %d,  t* = %d\n", Nn, input$d, input$T, input$tstar))
    cat(sprintf("ρ(W) = %.6f\n", stab$rhoW))
    cat(sprintf("α + β + γ·ρ(W) = %.6f  -> %s\n",
                stab$lhs, if (stab$ok) "OK (<1)" else "VIOLATES (<1 required)"))
  })
  output$stabilityInfo_diag <- renderPrint({
    stab <- stability()
    Nn <- N()
    cat(sprintf("N = %d  (d = %d),  T = %d,  t* = %d\n", Nn, input$d, input$T, input$tstar))
    cat(sprintf("ρ(W) = %.6f\n", stab$rhoW))
    cat(sprintf("α + β + γ·ρ(W) = %.6f  -> %s\n",
                stab$lhs, if (stab$ok) "OK (<1)" else "VIOLATES (<1 required)"))
    if (identical(input$omega_src, "upload_vec")) {
      if (is.null(input$omega_file)) cat("ω: upload mode (no file) → scalar fallback\n")
      else {
        vec <- omega_parsed()
        if (!is.null(vec)) cat(sprintf("ω JSON length = %d (current N=%d)\n", length(vec), Nn))
      }
    } else cat(sprintf("ω (scalar) = %.4f\n", input$omega_scalar))
    if (identical(input$W_src, "upload_W")) {
      if (is.null(input$W_file)) cat("W: upload mode (no file) → Toeplitz fallback\n")
      else {
        Wm <- W_parsed()
        if (!is.null(Wm)) cat(sprintf("W JSON size = %dx%d (current N=%d)%s\n",
                                      nrow(Wm), ncol(Wm), Nn,
                                      if (isTRUE(input$normalizeW_up)) " [row-normalized]" else ""))
      }
    } else {
      cat(sprintf("W: Toeplitz (diag=0), ρ=%.2f%s\n",
                  input$rhoW, if (isTRUE(input$normalizeW)) " [row-normalized]" else ""))
    }
    cat(sprintf("Animation: %s | FPS=%d | Loop=%s\n",
                if (isTRUE(playing())) "Playing" else "Paused",
                as.integer(input$fps),
                if (isTRUE(input$loop)) "Yes" else "No"))
    invisible()
  })
}

shinyApp(ui, server)
