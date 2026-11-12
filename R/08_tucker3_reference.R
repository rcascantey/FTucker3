# Extracted without numerical changes from M0 of the complete FTucker-3 thesis code.
# Original file: FTucker3_codigo_R_completo_DESCARGA(1).R
# These are library function definitions ONLY. Original experiment scripts are not run.

center_3way_tensor <- function(X) {
  d <- dim(X)
  I <- d[1]
  J <- d[2]
  K <- d[3]

  global <- mean(X)
  mi <- apply(X, 1, mean)
  mj <- apply(X, 2, mean)
  mk <- apply(X, 3, mean)

  Xc <- array(0, dim = d, dimnames = dimnames(X))

  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      for (k in seq_len(K)) {
        Xc[i, j, k] <- X[i, j, k] -
          mi[i] - mj[j] - mk[k] + 2 * global
      }
    }
  }

  list(
    X_centered = Xc,
    mean_global = global,
    mean_i = mi,
    mean_j = mj,
    mean_k = mk
  )
}

tucker3_single_start <- function(
    X,
    P,
    Q,
    R,
    max_iter = 300,
    tol = 1e-8,
    init = c("svd", "random")
) {
  init <- match.arg(init)

  d <- dim(X)

  if (init == "svd") {
    A <- svd(unfold_tensor(X, 1), nu = P, nv = 0)$u[, seq_len(P), drop = FALSE]
    B <- svd(unfold_tensor(X, 2), nu = Q, nv = 0)$u[, seq_len(Q), drop = FALSE]
    Cmat <- svd(unfold_tensor(X, 3), nu = R, nv = 0)$u[, seq_len(R), drop = FALSE]
  } else {
    A <- random_orthonormal(d[1], P)
    B <- random_orthonormal(d[2], Q)
    Cmat <- random_orthonormal(d[3], R)
  }

  prev <- Inf
  converged <- FALSE

  for (iter in seq_len(max_iter)) {

    ZA <- mode_n_product(X, t(B), 2)
    ZA <- mode_n_product(ZA, t(Cmat), 3)
    Anew <- svd(unfold_tensor(ZA, 1), nu = P, nv = 0)$u[, seq_len(P), drop = FALSE]

    ZB <- mode_n_product(X, t(Anew), 1)
    ZB <- mode_n_product(ZB, t(Cmat), 3)
    Bnew <- svd(unfold_tensor(ZB, 2), nu = Q, nv = 0)$u[, seq_len(Q), drop = FALSE]

    ZC <- mode_n_product(X, t(Anew), 1)
    ZC <- mode_n_product(ZC, t(Bnew), 2)
    Cnew <- svd(unfold_tensor(ZC, 3), nu = R, nv = 0)$u[, seq_len(R), drop = FALSE]

    G <- mode_n_product(X, t(Anew), 1)
    G <- mode_n_product(G, t(Bnew), 2)
    G <- mode_n_product(G, t(Cnew), 3)

    Xhat <- reconstruct_3way(G, Anew, Bnew, Cnew)

    err <- frobenius_norm(X - Xhat) / frobenius_norm(X)

    if (is.finite(prev) && abs(err - prev) < tol) {
      converged <- TRUE
      A <- Anew
      B <- Bnew
      Cmat <- Cnew
      break
    }

    A <- Anew
    B <- Bnew
    Cmat <- Cnew
    prev <- err
  }

  list(
    A = A,
    B = B,
    C = Cmat,
    G = G,
    X_hat = Xhat,
    relative_error = err,
    R2 = 1 - err^2,
    iterations = iter,
    converged = converged
  )
}

#' Fit classical three-way Tucker model
#'
#' Reference Tucker-3 alternating orthogonal estimation routine retained from the original controlled-experiment module.
#'
#' @param P Rank for the first structural mode.
#' @param Q Rank for the second structural mode.
#' @param R Rank for the third structural mode.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @param n_random_starts Additional random initializations.
#' @param X Three-way tensor to approximate.
#'
#' @return A list with loadings, core, reconstruction and errors.
#' @export
tucker3_fit <- function(
    X,
    P,
    Q,
    R,
    max_iter = 300,
    tol = 1e-8,
    n_random_starts = 5
) {
  fits <- list(
    tucker3_single_start(
      X, P, Q, R,
      max_iter = max_iter,
      tol = tol,
      init = "svd"
    )
  )

  if (n_random_starts > 0) {
    for (s in seq_len(n_random_starts)) {
      fits[[length(fits) + 1]] <- tucker3_single_start(
        X, P, Q, R,
        max_iter = max_iter,
        tol = tol,
        init = "random"
      )
    }
  }

  errors <- vapply(fits, function(z) z$relative_error, numeric(1))
  fits[[which.min(errors)]]
}

