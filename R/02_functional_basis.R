# Extracted without numerical changes from M0 of the complete FTucker-3 thesis code.
# Original file: FTucker3_codigo_R_completo_DESCARGA(1).R
# These are library function definitions ONLY. Original experiment scripts are not run.

#' Construct an orthonormal Fourier basis
#'
#' Evaluate the same orthonormal Fourier basis used in the original FTucker-3 thesis code.
#'
#' @param t Numeric evaluation times.
#' @param a Origin of the functional interval.
#' @param L Period of the Fourier basis.
#' @param N Number of Fourier basis functions (odd).
#'
#' @return A numeric matrix with one column per basis element.
#' @export
fourier_basis_matrix <- function(t, a = 0, L = 12, N = 5) {
  if (N < 1 || N %% 2 == 0) {
    stop("N must be odd: N = 2H + 1.")
  }

  H <- (N - 1) / 2
  Phi <- matrix(0, nrow = length(t), ncol = N)
  Phi[, 1] <- 1 / sqrt(L)

  if (H >= 1) {
    for (h in seq_len(H)) {
      Phi[, 2 * h] <- sqrt(2 / L) *
        cos(2 * pi * h * (t - a) / L)

      Phi[, 2 * h + 1] <- sqrt(2 / L) *
        sin(2 * pi * h * (t - a) / L)
    }
  }

  colnames(Phi) <- c(
    "Constant",
    as.vector(rbind(
      paste0("Cos", seq_len(H)),
      paste0("Sin", seq_len(H))
    ))
  )[seq_len(N)]

  Phi
}

weighted_ls_coef <- function(y, Phi, w = NULL) {
  if (is.null(w)) w <- rep(1, length(y))
  Wsqrt <- sqrt(w)

  Xw <- Phi * Wsqrt
  yw <- y * Wsqrt

  # qr.solve is numerically more stable than explicit matrix inversion.
  as.vector(qr.solve(Xw, yw))
}

reconstruct_curve <- function(coef, Phi) {
  as.vector(Phi %*% coef)
}

#' Select Fourier dimension by time-point cross-validation
#'
#' Implements the original leave-one-time-point-out cross-validation for selecting the number of Fourier basis functions.
#'
#' @param a Origin of the functional interval.
#' @param L Period of the Fourier basis.
#' @param X_std Four-way array (I, J, K, M) containing preprocessed functional observations.
#' @param t_grid Time grid corresponding to the fourth dimension.
#' @param candidates Candidate basis dimensions.
#' @param near_tie Near-tie tolerance for basis CV.
#'
#' @return A list containing selected dimension and a CV results table.
#' @export
select_fourier_dimension <- function(
    X_std,
    t_grid,
    candidates = c(3, 5, 7, 9),
    a = 0,
    L = 12,
    near_tie = 0.01
) {
  dims <- dim(X_std)
  I <- dims[1]
  J <- dims[2]
  K <- dims[3]
  M <- dims[4]

  candidates <- sort(unique(candidates))
  candidates <- candidates[
    candidates %% 2 == 1 &
      candidates >= 3 &
      candidates <= (M - 1)
  ]

  if (length(candidates) == 0) {
    stop("No admissible Fourier dimensions.")
  }

  results <- data.frame(
    N = candidates,
    CV_MSE = NA_real_,
    InSample_R2 = NA_real_
  )

  for (cc in seq_along(candidates)) {
    N <- candidates[cc]
    Phi_full <- fourier_basis_matrix(t_grid, a = a, L = L, N = N)

    cv_sse <- 0
    cv_n <- 0
    sse <- 0
    sst <- 0

    for (i in seq_len(I)) {
      for (j in seq_len(J)) {
        for (k in seq_len(K)) {

          y <- as.vector(X_std[i, j, k, ])

          # Full fit for descriptive R^2
          coef_full <- weighted_ls_coef(y, Phi_full)
          yhat_full <- reconstruct_curve(coef_full, Phi_full)

          sse <- sse + sum((y - yhat_full)^2)
          sst <- sst + sum(y^2)

          # Leave one time point out
          for (m in seq_len(M)) {
            keep <- setdiff(seq_len(M), m)

            coef_loo <- weighted_ls_coef(
              y = y[keep],
              Phi = Phi_full[keep, , drop = FALSE]
            )

            pred_m <- as.numeric(
              Phi_full[m, , drop = FALSE] %*% coef_loo
            )

            cv_sse <- cv_sse + (y[m] - pred_m)^2
            cv_n <- cv_n + 1
          }
        }
      }
    }

    results$CV_MSE[cc] <- cv_sse / cv_n
    results$InSample_R2[cc] <- 1 - sse / sst
  }

  min_cv <- min(results$CV_MSE)
  admissible <- results$N[
    results$CV_MSE <= min_cv * (1 + near_tie)
  ]

  N_selected <- min(admissible)

  list(
    selected = N_selected,
    table = results
  )
}

#' Expand functional three-way data into a coefficient tensor
#'
#' Projects each trajectory in a four-dimensional array (I, J, K, M) onto the Fourier basis, yielding a fourth-order coefficient tensor (I, J, K, N).
#'
#' @param a Origin of the functional interval.
#' @param L Period of the Fourier basis.
#' @param N Number of Fourier basis functions (odd).
#' @param X_std Four-way array (I, J, K, M) containing preprocessed functional observations.
#' @param t_grid Time grid corresponding to the fourth dimension.
#'
#' @return A list with elements C and Phi.
#' @export
build_coefficient_tensor <- function(
    X_std,
    t_grid,
    N,
    a = 0,
    L = 12
) {
  dims <- dim(X_std)
  I <- dims[1]
  J <- dims[2]
  K <- dims[3]

  Phi <- fourier_basis_matrix(t_grid, a = a, L = L, N = N)

  Ccoef <- array(
    0,
    dim = c(I, J, K, N),
    dimnames = c(
      dimnames(X_std)[1:3],
      list(basis = colnames(Phi))
    )
  )

  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      for (k in seq_len(K)) {
        Ccoef[i, j, k, ] <- weighted_ls_coef(
          as.vector(X_std[i, j, k, ]),
          Phi
        )
      }
    }
  }

  list(C = Ccoef, Phi = Phi)
}

