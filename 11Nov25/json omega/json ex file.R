library(jsonlite)

# Suppose d = 10, so N = 100
d <- 10
omega_vec <- rep(0.1, d * d)   # example ω vector of length N = d^2

# Write to JSON file (compact array form)
write_json(omega_vec, path = "omega_vector.json", pretty = FALSE, auto_unbox = TRUE)
