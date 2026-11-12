# Extracted without numerical changes from M0 of the complete FTucker-3 thesis code.
# Original file: FTucker3_codigo_R_completo_DESCARGA(1).R
# These are library function definitions ONLY. Original experiment scripts are not run.

#' Reconstruct curves from coefficients
#'
#' Evaluate the coefficient tensor over the rows of a supplied basis matrix.
#'
#' @param Ccoef Fourth-order coefficient tensor.
#' @param Phi Basis evaluation matrix.
#'
#' @return A fourth-order array of reconstructed functional observations.
#' @export
coefficients_to_functional_array <- function(Ccoef, Phi) {
  d <- dim(Ccoef)

  if (length(d) != 4) {
    stop("Ccoef must be a fourth-order tensor.")
  }

  I <- d[1]
  J <- d[2]
  K <- d[3]
  N <- d[4]
  M <- nrow(Phi)

  stopifnot(ncol(Phi) == N)

  # The reconstructed coefficient tensor may lose its dimnames during
  # successive mode-n products. Build a valid four-element dimnames list
  # explicitly, using existing labels when they are available.
  dn <- dimnames(Ccoef)

  get_labels <- function(idx, n) {
    if (!is.null(dn) &&
        length(dn) >= idx &&
        !is.null(dn[[idx]]) &&
        length(dn[[idx]]) == n) {
      return(dn[[idx]])
    }
    as.character(seq_len(n))
  }

  time_labels <- rownames(Phi)
  if (is.null(time_labels) || length(time_labels) != M) {
    time_labels <- as.character(seq_len(M))
  }

  output_dimnames <- list(
    mode1 = get_labels(1, I),
    mode2 = get_labels(2, J),
    mode3 = get_labels(3, K),
    time = time_labels
  )

  Xhat <- array(
    NA_real_,
    dim = c(I, J, K, M),
    dimnames = output_dimnames
  )

  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      for (k in seq_len(K)) {
        Xhat[i, j, k, ] <- as.vector(
          Phi %*% Ccoef[i, j, k, ]
        )
      }
    }
  }

  Xhat
}

#' Reconstruct the functional core
#'
#' Evaluate each functional-core coefficient series over a supplied basis matrix.
#'
#' @param Phi Basis evaluation matrix.
#' @param G Fourth-order core tensor.
#'
#' @return A fourth-order array with the functional core trajectories.
#' @export
functional_core_array <- function(G, Phi) {
  d <- dim(G)
  P <- d[1]
  Q <- d[2]
  R <- d[3]
  N <- d[4]
  M <- nrow(Phi)

  stopifnot(ncol(Phi) == N)

  H <- array(
    0,
    dim = c(P, Q, R, M)
  )

  for (p in seq_len(P)) {
    for (q in seq_len(Q)) {
      for (r in seq_len(R)) {
        H[p, q, r, ] <- as.vector(
          Phi %*% G[p, q, r, ]
        )
      }
    }
  }

  H
}

