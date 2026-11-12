# Extracted without numerical changes from M0 of the complete FTucker-3 thesis code.
# Original file: FTucker3_codigo_R_completo_DESCARGA(1).R
# These are library function definitions ONLY. Original experiment scripts are not run.

make_cell_folds <- function(I, J, K, Kfold = 5, seed = 1234) {
  set.seed(seed)
  n_cells <- I * J * K

  folds <- sample(rep(seq_len(Kfold), length.out = n_cells))
  array(folds, dim = c(I, J, K))
}

#' Select structural ranks by cell-wise cross-validation
#'
#' Select P, Q, and R using training-cell centering, withheld whole functional cells, and the one-standard-error parsimony rule.
#'
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @param C_raw Raw coefficient tensor.
#' @param P_candidates Candidate ranks for mode 1.
#' @param Q_candidates Candidate ranks for mode 2.
#' @param R_candidates Candidate ranks for mode 3.
#' @param Kfold Number of cellwise CV folds.
#' @param seed Random seed used to generate folds.
#'
#' @return A list with the selected ranks, CV table, folds, and threshold.
#' @export
rank_cv_ftucker <- function(
    C_raw,
    P_candidates,
    Q_candidates,
    R_candidates,
    Kfold = 5,
    max_iter = 150,
    tol = 1e-7,
    seed = 1234
) {
  d <- dim(C_raw)
  I <- d[1]
  J <- d[2]
  K <- d[3]
  N <- d[4]

  fold_array <- make_cell_folds(I, J, K, Kfold, seed)

  grid <- expand.grid(
    P = P_candidates,
    Q = Q_candidates,
    R = R_candidates
  ) %>%
    arrange(P, Q, R)

  out <- vector("list", nrow(grid))

  for (g in seq_len(nrow(grid))) {
    P <- grid$P[g]
    Q <- grid$Q[g]
    R <- grid$R[g]

    fold_mse <- numeric(Kfold)

    cat(
      sprintf(
        "Rank CV: P=%d, Q=%d, R=%d (%d/%d)\n",
        P, Q, R, g, nrow(grid)
      )
    )

    for (f in seq_len(Kfold)) {
      mask_train <- array(1, dim = c(I, J, K))
      mask_train[fold_array == f] <- 0

      # Training-only marginal centering
      cent <- center_coefficient_tensor(
        Ccoef = C_raw,
        mask3 = mask_train
      )

      Cc <- cent$C_centered

      fit <- ftuckals3(
        Ccenter = Cc,
        P = P, Q = Q, R = R,
        mask3 = mask_train,
        max_iter = max_iter,
        tol = tol,
        n_random_starts = 0,
        verbose = FALSE
      )

      mask_val4 <- array(
        rep((1 - mask_train), times = N),
        dim = dim(C_raw)
      )

      resid_val <- Cc - fit$C_hat
      fold_mse[f] <- mean(
        resid_val[mask_val4 == 1]^2,
        na.rm = TRUE
      )
    }

    model_df <-
      P * (I - P) +
      Q * (J - Q) +
      R * (K - R) +
      P * Q * R * N

    out[[g]] <- data.frame(
      P = P,
      Q = Q,
      R = R,
      CoreSize = P * Q * R,
      ModelDF = model_df,
      CV_MSE = mean(fold_mse),
      CV_SD = stats::sd(fold_mse),
      CV_SE = stats::sd(fold_mse) / sqrt(Kfold)
    )
  }

  results <- bind_rows(out)

  # One-standard-error parsimonious rule
  imin <- which.min(results$CV_MSE)
  threshold <- results$CV_MSE[imin] + results$CV_SE[imin]

  eligible <- results %>%
    filter(CV_MSE <= threshold) %>%
    arrange(ModelDF, CV_MSE, P, Q, R)

  selected <- eligible[1, , drop = FALSE]

  list(
    selected = selected,
    table = results,
    folds = fold_array,
    threshold = threshold
  )
}

