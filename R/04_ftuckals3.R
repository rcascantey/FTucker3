# Extracted without numerical changes from M0 of the complete FTucker-3 thesis code.
# Original file: FTucker3_codigo_R_completo_DESCARGA(1).R
# These are library function definitions ONLY. Original experiment scripts are not run.

ftuckals3_single_start <- function(
    Ccenter,
    P,
    Q,
    R,
    mask3 = NULL,
    max_iter = 300,
    tol = 1e-8,
    init = c("svd", "random"),
    verbose = FALSE
) {
  init <- match.arg(init)

  d <- dim(Ccenter)
  I <- d[1]
  J <- d[2]
  K <- d[3]
  N <- d[4]

  if (is.null(mask3)) mask3 <- array(1, dim = c(I, J, K))
  mask4 <- array(rep(mask3, times = N), dim = d)

  Cobs <- Ccenter
  Cobs[mask4 == 0] <- NA_real_

  Cfilled <- Cobs
  Cfilled[is.na(Cfilled)] <- 0

  if (init == "svd") {
    A <- svd(unfold_tensor(Cfilled, 1), nu = P, nv = 0)$u[, seq_len(P), drop = FALSE]
    B <- svd(unfold_tensor(Cfilled, 2), nu = Q, nv = 0)$u[, seq_len(Q), drop = FALSE]
    Cmat <- svd(unfold_tensor(Cfilled, 3), nu = R, nv = 0)$u[, seq_len(R), drop = FALSE]
  } else {
    A <- random_orthonormal(I, P)
    B <- random_orthonormal(J, Q)
    Cmat <- random_orthonormal(K, R)
  }

  previous_error <- Inf
  converged <- FALSE

  for (iter in seq_len(max_iter)) {

    # Current reconstruction used to impute withheld cells
    Gcur <- mode_n_product(Cfilled, t(A), 1)
    Gcur <- mode_n_product(Gcur, t(B), 2)
    Gcur <- mode_n_product(Gcur, t(Cmat), 3)

    Cpred_cur <- reconstruct_4way(Gcur, A, B, Cmat)

    Cwork <- Cobs
    Cwork[mask4 == 0] <- Cpred_cur[mask4 == 0]

    # Update A
    ZA <- mode_n_product(Cwork, t(B), 2)
    ZA <- mode_n_product(ZA, t(Cmat), 3)
    Anew <- svd(
      unfold_tensor(ZA, 1),
      nu = P,
      nv = 0
    )$u[, seq_len(P), drop = FALSE]

    # Update B
    ZB <- mode_n_product(Cwork, t(Anew), 1)
    ZB <- mode_n_product(ZB, t(Cmat), 3)
    Bnew <- svd(
      unfold_tensor(ZB, 2),
      nu = Q,
      nv = 0
    )$u[, seq_len(Q), drop = FALSE]

    # Update C
    ZC <- mode_n_product(Cwork, t(Anew), 1)
    ZC <- mode_n_product(ZC, t(Bnew), 2)
    Cnew <- svd(
      unfold_tensor(ZC, 3),
      nu = R,
      nv = 0
    )$u[, seq_len(R), drop = FALSE]

    # Core and reconstruction
    Gnew <- mode_n_product(Cwork, t(Anew), 1)
    Gnew <- mode_n_product(Gnew, t(Bnew), 2)
    Gnew <- mode_n_product(Gnew, t(Cnew), 3)

    Chat <- reconstruct_4way(Gnew, Anew, Bnew, Cnew)

    rel_error <- frobenius_norm(Ccenter - Chat, mask4) /
      frobenius_norm(Ccenter, mask4)

    change <- abs(rel_error - previous_error)

    if (verbose && (iter == 1 || iter %% 10 == 0)) {
      cat(
        sprintf(
          "iter=%d | observed relative error=%.10f | change=%.3e\n",
          iter, rel_error, change
        )
      )
    }

    A <- Anew
    B <- Bnew
    Cmat <- Cnew
    Cfilled <- Cwork

    if (is.finite(previous_error) && change < tol) {
      converged <- TRUE
      break
    }

    previous_error <- rel_error
  }

  list(
    A = A,
    B = B,
    C = Cmat,
    G = Gnew,
    C_hat = Chat,
    relative_error = rel_error,
    R2_coeff = 1 - rel_error^2,
    iterations = iter,
    converged = converged,
    init = init
  )
}

#' Estimate FTucker-3 by alternating orthogonal updates
#'
#' Fit the original FTUCKALS-3 algorithm on a four-way coefficient tensor, reducing the first three modes while preserving the functional-basis mode.
#'
#' @param mask3 Optional structural-cell observation mask.
#' @param Ccenter Centered coefficient tensor.
#' @param P Rank for the first structural mode.
#' @param Q Rank for the second structural mode.
#' @param R Rank for the third structural mode.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @param n_random_starts Additional random initializations.
#' @param verbose Print progress when TRUE.
#'
#' @return A list with A, B, C, functional core G, coefficient reconstruction C_hat, and convergence measures.
#' @export
ftuckals3 <- function(
    Ccenter,
    P,
    Q,
    R,
    mask3 = NULL,
    max_iter = 300,
    tol = 1e-8,
    n_random_starts = 0,
    verbose = FALSE
) {
  fits <- list()

  # Rational/SVD initialization
  fits[[1]] <- ftuckals3_single_start(
    Ccenter = Ccenter,
    P = P, Q = Q, R = R,
    mask3 = mask3,
    max_iter = max_iter,
    tol = tol,
    init = "svd",
    verbose = verbose
  )

  if (n_random_starts > 0) {
    for (s in seq_len(n_random_starts)) {
      fits[[length(fits) + 1]] <- ftuckals3_single_start(
        Ccenter = Ccenter,
        P = P, Q = Q, R = R,
        mask3 = mask3,
        max_iter = max_iter,
        tol = tol,
        init = "random",
        verbose = FALSE
      )
    }
  }

  errors <- vapply(fits, function(z) z$relative_error, numeric(1))
  best <- which.min(errors)

  out <- fits[[best]]
  out$all_start_errors <- errors
  out$best_start <- best
  out
}

