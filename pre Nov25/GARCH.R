library(spGarch)

# 1st example (spatial GARCH)
##########################

# parameters
rho    <- 0.5
lambda <- 0.3
alpha  <- 1
d      <- 5

nblist <- cell2nb(d, d, type = "rook") # lattice process with Rook's contiguity matrix
M_1 <- nb2mat(nblist)
M_2 <- M_1

# simulation
Y <- sim.spGARCH(rho = rho, lambda = lambda, alpha = alpha,
                 W1 = M_1, W2 = M_2, type = "spGARCH")

# visualization
image(1:d, 1:d, array(Y, dim = c(d, d)),
      xlab = expression(s[1]), ylab = expression(s[2]))


# 2nd example (exponential spatial GARCH)
##########################

# parameters
rho    <- 0.5
lambda <- 0.3
alpha  <- 1
zeta   <- 0.5
theta  <- 0.5
d      <- 5

nblist <- cell2nb(d, d, type = "rook") # lattice process with Rook's contiguity matrix
M_1 <- nb2mat(nblist)
M_2 <- M_1

# simulation
Y <- sim.spGARCH(rho = rho, lambda = lambda, alpha = alpha,
                 W1 = M_1, W2 = M_2, zeta = zeta, theta = theta,
                 type = "e-spGARCH")

# visualization
image(1:d, 1:d, array(Y, dim = c(d, d)),
      xlab = expression(s[1]), ylab = expression(s[2]))



# Real data example 
##########################

