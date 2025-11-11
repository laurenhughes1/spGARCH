#-----------------------------------------------
# Purely spatial GARCH/ARCH model
#-----------------------------------------------

# Create Toeplitz weight matrix using base::toeplitz()
make_toeplitz <- function(N, rho = 0.5, normalize = TRUE) {
  # first column for toeplitz (decaying with distance)
  first_col <- rho^(0:(N - 1))
  W <- toeplitz(first_col)
  if (normalize) {
    W <- W / rowSums(W)
  }
  return(W)
}

# Spatial variance equation:
# (I - λW) σ² = ω + α p + γ W p
spatial_garch_once <- function(eps, omega, alpha, gamma, W, lambda = 0) {
  N <- length(eps)
  p <- eps^2
  rhs <- omega + alpha * p + gamma * (W %*% p)
  if (lambda == 0) {
    sigma2 <- rhs
  } else {
    sigma2 <- solve(diag(N) - lambda * W, rhs)
  }
  as.vector(sigma2)
}

#-----------------------------------------------
# Example (deterministic)
#-----------------------------------------------
set.seed(123)

N <- 25                      # number of spatial units (perfect square)
omega <- rep(0.1, N)
alpha <- 0.2
gamma <- 0.1
lambda <- 0.3                # 0 for spatial-ARCH; >0 for contemporaneous feedback
eps_fixed <- rnorm(N)        # pre-simulated white noise

# Build Toeplitz spatial weights
W <- make_toeplitz(N, rho = 0.5, normalize = TRUE)

# Compute one-shot spatial conditional variances
sigma2 <- spatial_garch_once(eps_fixed, omega, alpha, gamma, W, lambda)

#-----------------------------------------------
# Visualization
#-----------------------------------------------

# reshape into square grid
n_side <- sqrt(N)
Z <- matrix(sigma2, nrow = n_side, byrow = TRUE)

# Heatmap
image(
  1:n_side, 1:n_side, Z,
  col = heat.colors(40),
  main = expression("Heatmap of " * sigma^2 * " (purely spatial)"),
  xlab = "Spatial X", ylab = "Spatial Y"
)

# Contour plot
filled.contour(
  1:n_side, 1:n_side, Z,
  color.palette = terrain.colors,
  main = expression("Contour of " * sigma^2 * " (purely spatial)"),
  xlab = "Spatial X", ylab = "Spatial Y"
)
