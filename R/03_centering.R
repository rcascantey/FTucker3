# Extracted without numerical changes from M0 of the complete FTucker-3 thesis code.
# Original file: FTucker3_codigo_R_completo_DESCARGA(1).R
# These are library function definitions ONLY. Original experiment scripts are not run.

#' Center the coefficient tensor
#'
#' Apply the marginal additive centering specified in the original FTucker-3 implementation.
#'
#' @param Ccoef Fourth-order coefficient tensor.
#' @param mask3 Optional structural-cell observation mask.
#'
#' @return A list with centered coefficients and marginal means.
#' @export
center_coefficient_tensor <- function(Ccoef, mask3 = NULL) {
  d <- dim(Ccoef)
  I <- d[1]
  J <- d[2]
  K <- d[3]
  N <- d[4]

  if (is.null(mask3)) {
    mask3 <- array(1, dim = c(I, J, K))
  }

  Ccenter <- array(NA_real_, dim = d, dimnames = dimnames(Ccoef))

  mean_global <- numeric(N)
  mean_i <- matrix(NA_real_, I, N)
  mean_j <- matrix(NA_real_, J, N)
  mean_k <- matrix(NA_real_, K, N)

  for (n in seq_len(N)) {
    Xn <- Ccoef[, , , n]

    obs_global <- Xn[mask3 == 1]
    mean_global[n] <- mean(obs_global, na.rm = TRUE)

    for (i in seq_len(I)) {
      tmp <- Xn[i, , ]
      mm <- mask3[i, , ]
      mean_i[i, n] <- mean(tmp[mm == 1], na.rm = TRUE)
    }

    for (j in seq_len(J)) {
      tmp <- Xn[, j, ]
      mm <- mask3[, j, ]
      mean_j[j, n] <- mean(tmp[mm == 1], na.rm = TRUE)
    }

    for (k in seq_len(K)) {
      tmp <- Xn[, , k]
      mm <- mask3[, , k]
      mean_k[k, n] <- mean(tmp[mm == 1], na.rm = TRUE)
    }

    if (anyNA(mean_i[, n]) || anyNA(mean_j[, n]) || anyNA(mean_k[, n])) {
      stop(
        "A CV fold removed all observations from at least one marginal level. ",
        "Use fewer folds or a different fold allocation."
      )
    }

    for (i in seq_len(I)) {
      for (j in seq_len(J)) {
        for (k in seq_len(K)) {
          Ccenter[i, j, k, n] <-
            Ccoef[i, j, k, n] -
            mean_i[i, n] -
            mean_j[j, n] -
            mean_k[k, n] +
            2 * mean_global[n]
        }
      }
    }
  }

  list(
    C_centered = Ccenter,
    mean_global = mean_global,
    mean_i = mean_i,
    mean_j = mean_j,
    mean_k = mean_k
  )
}

#' Restore additive marginal effects
#'
#' Reverse the additive centering transformation applied to the coefficient tensor.
#'
#' @param Ccenter Centered coefficient tensor.
#' @param centering Marginal centering parameters returned by center_coefficient_tensor().
#'
#' @return The coefficient tensor after reversing marginal centering.
#' @export
decenter_coefficient_tensor <- function(Ccenter, centering) {
  d <- dim(Ccenter)
  I <- d[1]
  J <- d[2]
  K <- d[3]
  N <- d[4]

  Cout <- array(NA_real_, dim = d, dimnames = dimnames(Ccenter))

  for (n in seq_len(N)) {
    for (i in seq_len(I)) {
      for (j in seq_len(J)) {
        for (k in seq_len(K)) {
          Cout[i, j, k, n] <-
            Ccenter[i, j, k, n] +
            centering$mean_i[i, n] +
            centering$mean_j[j, n] +
            centering$mean_k[k, n] -
            2 * centering$mean_global[n]
        }
      }
    }
  }

  Cout
}

