#' The SOPS procedure
#'
#' @param X design points, should be standardized in (0, 1)
#' @param Y responses
#' @param h_seq candidates for the bandwidth
#' @param n_ss times of sample-splitting for bandwidth selection
#' @param kappa controlling the size of removed neighborhood of detected jumps
#' @param q nominal FDR level
#' @param thres.1 whether to increment the number of negatives by 1 to determine the threshold
#'
#' @return a list contains the number, locations and sizes of detected jumps
#' @export
#'
#' @examples NA
SOPS <- function(X, Y, h_seq = NULL, n_ss = 20, kappa = 1, q = 0.2, thres.1 = TRUE) {

  # reordering ----
  index_X <- sort(X, decreasing = FALSE, index.return = TRUE)$ix
  X <- X[index_X]
  Y <- Y[index_X]
  x_seq <- X

  n <- length(Y)

  # select h ----
  if (length(h_seq) == 1) {
    h_opt <- h_seq
  } else {
    res_COPS <- COPS(X, Y, h_seq, n_ss, kappa)
    h_opt <- res_COPS$h_opt
  }

  # SOPS ----

  # ++ OPS ----
  tt <- 1:n
  index_O <- tt[1:n %% 2 == 1]
  index_E <- tt[1:n %% 2 == 0]
  X_O <- X[index_O]
  X_E <- X[index_E]
  Y_O <- Y[index_O]
  Y_E <- Y[index_E]

  h <- h_opt
  L <- ceiling(1 / (kappa * h))  # (1-2h)/(kappa * h)

  # ++ O -> E ----
  res_O <- val(x_seq, Y_O, X_O, Y_E, X_E, h, kappa, L)

  # ++ W ----
  W <- n * h * res_O$dh_O * res_O$dh_E
  sFDP_temp <- sFDP(W, q)
  if (thres.1) {
    FDP <- sFDP_temp$FDP[, 2]
  } else {
    FDP <- sFDP_temp$FDP[, 1]
  }
  ok <- which(FDP <= q)
  if (length(ok) > 0) {
    threshold <- sFDP_temp$s[ok[1]]
  } else {
    threshold <- Inf
  }
  P_t <- W >= threshold

  list(
    number = sum(P_t),
    locations = res_O$th_O[P_t],
    sizes = res_O$dh_O[P_t],
    h_opt = h_opt
  )
}
