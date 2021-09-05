#' The COPS procedure
#'
#' @param X design points, should be standardized in (0, 1)
#' @param Y responses
#' @param h_seq candidates for the bandwidth
#' @param n_ss times of sample-splitting for bandwidth selection
#' @param kappa controlling the size of removed neighborhood of detected jumps
#'
#' @return a list contains the number, locations and sizes of detected jumps
#' @export
#'
#' @examples NA
COPS <- function(X, Y, h_seq = NULL, n_ss = 20, kappa = 1) {

  # reordering ----
  index_X <- sort(X, decreasing = FALSE, index.return = TRUE)$ix
  X <- X[index_X]
  Y <- Y[index_X]
  x_seq <- X

  n <- length(Y)

  if (length(h_seq) == 1) {

    res <- COPS_h(X, Y, h_seq, kappa)
    return(res)

  } else {

    # h_seq ----
    if (is.null(h_seq)) {
      h_min <- 0.01
      h_max <- 0.4
      a <- 0.8
      h_seq <- h_max * a^(0:floor(log(h_min / h_max, a)))
    } else {
      h_seq <- sort(h_seq, decreasing = TRUE)
    }
    n_h_seq <- length(h_seq)

    # select h ----
    jumps_n_ss <- matrix(NA, n_h_seq, n_ss)
    # jumps_loc_ss <- matrix(list(), n_h_seq, n_ss)
    te_H <- matrix(NA, n_h_seq, n_ss)
    for (i_ss in 1:n_ss) {

      # ++ OPS ----
      IDX <- OPS(n, prop = c(1, 1, 1))
      index_O <- IDX[1, IDX[1, ] > 0]
      index_E <- IDX[2, IDX[2, ] > 0]
      index_H <- IDX[3, IDX[3, ] > 0]
      X_O <- X[index_O]
      X_E <- X[index_E]
      X_H <- X[index_H]
      Y_O <- Y[index_O]
      Y_E <- Y[index_E]
      Y_H <- Y[index_H]

      for (i_h in 1:n_h_seq) {

        h <- h_seq[i_h]
        L <- ceiling(1 / (kappa * h))  # (1-2h)/(kappa * h)

        # ++ COPS-O ----
        res_O <- val(x_seq, Y_O, X_O, Y_E, X_E, h, kappa, L)
        te_E <- res_O$val_E
        jumps_n_O <- which.min(te_E)
        jumps_n_ss[i_h, i_ss] <- jumps_n_O
        jumps_loc_O <- as.vector(res_O$th_O[1:jumps_n_O])
        # jumps_loc_ss[i_h, i_ss][[1]] <- jumps_loc_O
        jumps_loc_O_sort <- sort(jumps_loc_O, index.return = T)
        th_O <- jumps_loc_O_sort$x
        dh_O <- res_O$dh_O[jumps_loc_O_sort$ix]
        res_H <- RSS_h(th_O, dh_O, Y_O, X_O, Y_H, X_H, h)
        te_H[i_h, i_ss] <- res_H$RSS
      }
    }
    te_H_h <- apply(te_H, 2, which.min)
    idx_h <- min(te_H_h)
    h_opt <- h_seq[idx_h]

    res <- COPS_h(X, Y, h_opt, kappa)
    res$h_opt <- h_opt
    return(res)
  }
}
