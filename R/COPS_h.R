COPS_h <- function(X, Y, h, kappa = 1) {

  # reordering ----
  index_X <- sort(X, decreasing = FALSE, index.return = TRUE)$ix
  X <- X[index_X]
  Y <- Y[index_X]
  x_seq <- X

  n <- length(Y)

  L <- ceiling(1 / (kappa * h))  # (1-2h)/(kappa * h)

  # OPS ----
  tt <- 1:n
  index_O <- tt[1:n %% 2 == 1]
  index_E <- tt[1:n %% 2 == 0]
  X_O <- X[index_O]
  X_E <- X[index_E]
  Y_O <- Y[index_O]
  Y_E <- Y[index_E]

  # COPS-O ----
  res_O <- val(x_seq, Y_O, X_O, Y_E, X_E, h, kappa, L)
  te_E <- res_O$val_E
  # COPS-E ----
  res_E <- val(x_seq, Y_E, X_E, Y_O, X_O, h, kappa, L)
  te_O <- res_E$val_E
  # COPS-CV ----
  n_OE <- min(length(te_E), length(te_O))
  te_OE <- te_E[1:n_OE] + te_O[1:n_OE]
  jumps_n_est <- which.min(te_OE)

  # refit ----
  res <- detect(x_seq, Y, X, h, kappa, jumps_n_est)
  locations <- res$jump_loc
  sizes <- res$jump_size

  list(
    number = jumps_n_est,
    locations = locations,
    sizes = sizes
  )
}
