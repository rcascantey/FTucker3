# Extracted without numerical changes from M0 of the complete FTucker-3 thesis code.
# Original file: FTucker3_codigo_R_completo_DESCARGA(1).R
# These are library function definitions ONLY. Original experiment scripts are not run.

frobenius_norm <- function(X, mask = NULL) {
  if (is.null(mask)) {
    return(sqrt(sum(X^2, na.rm = TRUE)))
  }
  sqrt(sum(X[mask == 1]^2, na.rm = TRUE))
}

unfold_tensor <- function(X, mode) {
  d <- dim(X)
  stopifnot(mode >= 1, mode <= length(d))
  perm <- c(mode, setdiff(seq_along(d), mode))
  Xp <- aperm(X, perm)
  matrix(Xp, nrow = d[mode])
}

mode_n_product <- function(X, M, mode) {
  # X: tensor
  # M: matrix of dimension r x dim(X)[mode]
  d <- dim(X)
  stopifnot(ncol(M) == d[mode])

  perm <- c(mode, setdiff(seq_along(d), mode))
  Xp <- aperm(X, perm)
  Xmat <- matrix(Xp, nrow = d[mode])

  Ymat <- M %*% Xmat
  d_new_perm <- c(nrow(M), d[setdiff(seq_along(d), mode)])
  Yp <- array(Ymat, dim = d_new_perm)

  aperm(Yp, order(perm))
}

reconstruct_4way <- function(G, A, B, Cmat) {
  Xhat <- mode_n_product(G, A, 1)
  Xhat <- mode_n_product(Xhat, B, 2)
  Xhat <- mode_n_product(Xhat, Cmat, 3)
  Xhat
}

reconstruct_3way <- function(G, A, B, Cmat) {
  Xhat <- mode_n_product(G, A, 1)
  Xhat <- mode_n_product(Xhat, B, 2)
  Xhat <- mode_n_product(Xhat, Cmat, 3)
  Xhat
}

random_orthonormal <- function(n, r) {
  if (r > n) stop("r cannot exceed n.")
  Q <- qr.Q(qr(matrix(stats::rnorm(n * r), nrow = n, ncol = r)))
  Q[, seq_len(r), drop = FALSE]
}

random_centered_orthonormal <- function(n, r) {
  # Orthonormal columns constrained to be orthogonal to the constant vector.
  if (r >= n) stop("For centered orthonormal factors, r must be < n.")
  u0 <- rep(1 / sqrt(n), n)
  Z <- matrix(stats::rnorm(n * n), nrow = n)
  Z <- Z - u0 %*% (t(u0) %*% Z)
  Q <- qr.Q(qr(Z))
  Q[, seq_len(r), drop = FALSE]
}

#' Rotation-invariant subspace similarity
#'
#' Calculate normalized squared Frobenius norm of cross-loadings for equal-dimensional orthonormal loading spaces.
#'
#' @param U_true First orthonormal loading matrix.
#' @param U_est Second orthonormal loading matrix.
#'
#' @return One numeric similarity value.
#' @export
subspace_similarity <- function(U_true, U_est) {
  r <- ncol(U_true)
  if (ncol(U_est) != r) stop("Subspaces must have equal dimensions.")
  sum((t(U_true) %*% U_est)^2) / r
}

