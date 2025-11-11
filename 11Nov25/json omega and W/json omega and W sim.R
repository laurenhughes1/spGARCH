# ============================================================
# Generate and save JSON: omega (length N=d^2) and W (N x N)
# - W: geometric-decay Toeplitz, diag=0, row-normalized
# - omega: smooth spatial variation, positive and bounded
# Files written: "omega_vector.json", "W_matrix.json"
# ============================================================

library(jsonlite)

set.seed(123)

# ---- User parameters ----
d      <- 10              # grid size; N = d^2
N      <- d * d
rhoW   <- 0.6             # decay parameter for geometric Toeplitz (0 <= rho < 1)
omega_min <- 0.05         # lower bound for omega
omega_max <- 0.15         # upper bound for omega
pretty_json <- TRUE       # pretty-print JSON for readability

# ---- Helper: Toeplitz W with diag=0 and row-normalization ----
make_toeplitz_W <- function(N, rho, normalize = TRUE) {
  stopifnot(N >= 2, rho >= 0, rho < 1)
  first_col <- rho^(0:(N - 1))
  W <- toeplitz(first_col)
  diag(W) <- 0
  if (normalize) {
    rs <- rowSums(W)
    rs[rs == 0] <- 1
    W <- W / rs
  }
  W
}

# ---- Build W ----
W <- make_toeplitz_W(N = N, rho = rhoW, normalize = TRUE)

# ---- Build a spatially smooth omega of length N = d^2 ----
# Create a gentle gradient over the 2D grid, then bound to [omega_min, omega_max]
coords <- expand.grid(y = seq_len(d), x = seq_len(d))
x_scaled <- (coords$x - 1) / (d - 1)
y_scaled <- (coords$y - 1) / (d - 1)

# Smooth field: center + gradients + tiny noise
center_val <- (omega_min + omega_max) / 2
amp        <- (omega_max - omega_min) / 2 * 0.9  # leave margin to stay within bounds
omega_raw  <- center_val +
  amp * (0.6 * x_scaled + 0.4 * y_scaled) +    # gentle slope across grid
  rnorm(N, sd = 0.005)                          # tiny noise

# Clip to [omega_min, omega_max]
omega <- pmin(omega_max, pmax(omega_min, omega_raw))

# ---- (Optional) quick sanity checks ----
cat(sprintf("N = %d (d = %d)\n", N, d))
cat(sprintf("omega range: [%.4f, %.4f]\n", min(omega), max(omega)))
cat(sprintf("W dims: %d x %d | row sums (min, max): (%.4f, %.4f)\n",
            nrow(W), ncol(W), min(rowSums(W)), max(rowSums(W))))

# ---- Save to JSON ----
write_json(omega, path = "omega_vector.json", pretty = pretty_json, auto_unbox = TRUE)
write_json(W,     path = "W_matrix.json",     pretty = pretty_json, auto_unbox = TRUE)

cat('Wrote files:\n  - "omega_vector.json"\n  - "W_matrix.json"\n')
