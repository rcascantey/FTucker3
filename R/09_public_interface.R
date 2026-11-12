#' Fit the FTucker-3 model to three-way functional data
#'
#' Convenience interface for the unmodified Fourier expansion, coefficient
#' centering, and FTUCKALS-3 functions extracted from the original doctoral
#' implementation. Functional observations are supplied as an array with
#' dimensions (I, J, K, M). If variables have different scales, standardize
#' the data before calling this function.
#'
#' @param X Numeric four-dimensional array with dimensions (I,J,K,M).
#' @param t_grid Numeric vector of M observation times.
#' @param ranks Three positive integers (P,Q,R), each not exceeding the
#'   respective structural mode dimension.
#' @param N Odd positive number of Fourier basis functions; must not exceed M.
#' @param L Period of the Fourier basis.
#' @param a Origin of the Fourier interval.
#' @param center Whether to apply marginal coefficient centering.
#' @param max_iter Maximum number of FTUCKALS-3 iterations.
#' @param tol Convergence tolerance for the iterative algorithm.
#' @param n_random_starts Additional random initializations.
#' @param verbose Print iterative progress if TRUE.
#' @return An FTucker3Fit object containing the coefficient tensor, loadings,
#'   functional core, reconstructed coefficient tensor and reconstructed data.
#' @examples
#' set.seed(23)
#' X <- array(rnorm(3 * 3 * 3 * 8), c(3, 3, 3, 8))
#' fit <- ftucker3(X, t_grid = 0:7, ranks = c(2, 2, 1), N = 3, L = 8,
#'                  max_iter = 10)
#' dim(fit$reconstructed)
#' @export
ftucker3 <- function(X, t_grid, ranks, N = 3, L, a = 0,
                     center = TRUE, max_iter = 300, tol = 1e-8,
                     n_random_starts = 0, verbose = FALSE) {
  d <- dim(X)
  if (!is.numeric(X) || length(d) != 4L || any(d < 2L) ||
      any(!is.finite(X))) {
    stop("X must be a finite, numeric 4-way array, with mode sizes >= 2.")
  }
  if (length(t_grid) != d[4] || any(!is.finite(t_grid))) {
    stop("t_grid must contain one finite value per functional time point.")
  }
  if (missing(L) || length(L) != 1L || !is.finite(L) || L <= 0) {
    stop("Provide a finite positive period L (e.g., L = 12 for monthly data).")
  }
  if (length(N) != 1L || !is.finite(N) || N != as.integer(N) ||
      N < 1L || N %% 2L == 0L || N > d[4]) {
    stop("N must be a positive odd integer not exceeding the number of times.")
  }
  if (length(ranks) != 3L || any(!is.finite(ranks)) ||
      any(ranks != as.integer(ranks)) || any(ranks < 1L) ||
      any(ranks > d[1:3])) {
    stop("ranks must contain three integers between 1 and mode dimensions.")
  }
  if (length(center) != 1L || is.na(center) || !is.logical(center)) {
    stop("center must be TRUE or FALSE.")
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
      length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter != as.integer(max_iter) || max_iter < 1L ||
      length(n_random_starts) != 1L || !is.finite(n_random_starts) ||
      n_random_starts != as.integer(n_random_starts) || n_random_starts < 0L) {
    stop("Invalid optimization options.")
  }
  built <- build_coefficient_tensor(X, t_grid = t_grid, N = N, a = a, L = L)
  centering <- if (center) center_coefficient_tensor(built$C) else NULL
  Ccenter <- if (center) centering$C_centered else built$C
  fit <- ftuckals3(Ccenter, P = ranks[1], Q = ranks[2], R = ranks[3],
                   max_iter = max_iter, tol = tol,
                   n_random_starts = n_random_starts, verbose = verbose)
  C_hat <- if (center) decenter_coefficient_tensor(fit$C_hat, centering) else fit$C_hat
  out <- list(
    loadings = list(A = fit$A, B = fit$B, C = fit$C),
    core = fit$G,
    functional_core = functional_core_array(fit$G, built$Phi),
    coefficient_tensor = built$C,
    reconstructed_coefficients = C_hat,
    reconstructed = coefficients_to_functional_array(C_hat, built$Phi),
    basis = built$Phi,
    centering = centering,
    optimization = fit,
    time_grid = t_grid,
    ranks = as.integer(ranks), N = N, L = L, a = a
  )
  class(out) <- "FTucker3Fit"
  out
}

#' Predict functional trajectories from an FTucker3Fit object
#'
#' @param object An FTucker3Fit object created with ftucker3().
#' @param t_grid Optional new evaluation times, defaulting to original times.
#' @param ... Reserved for future use.
#' @return A four-dimensional array of reconstructed functions.
#' @export
predict.FTucker3Fit <- function(object, t_grid = object$time_grid, ...) {
  Phi <- fourier_basis_matrix(t_grid, a = object$a, L = object$L, N = object$N)
  coefficients_to_functional_array(object$reconstructed_coefficients, Phi)
}
