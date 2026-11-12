# Extracted without numerical changes from M0 of the complete FTucker-3 thesis code.
# Original file: FTucker3_codigo_R_completo_DESCARGA(1).R
# These are library function definitions ONLY. Original experiment scripts are not run.

#' Compute structural-mode singular spectra
#'
#' Compute singular-value explained-variation proportions for one mode of the coefficient tensor.
#'
#' @param Ccenter Centered coefficient tensor.
#' @param mode Tensor mode index (1, 2 or 3).
#'
#' @return A data frame of explained and cumulative percentages.
#' @export
mode_spectrum <- function(Ccenter, mode) {
  s2 <- svd(unfold_tensor(Ccenter, mode), nu = 0, nv = 0)$d^2
  data.frame(
    Component = seq_along(s2),
    Explained = 100 * s2 / sum(s2),
    Cumulative = 100 * cumsum(s2) / sum(s2)
  )
}

#' Summarize energy of functional core interactions
#'
#' Compute squared-coefficient energy for each combination of structural core indices.
#'
#' @param G Fourth-order core tensor.
#'
#' @return A data frame of core interaction energies.
#' @export
core_energy_table <- function(G) {
  d <- dim(G)
  P <- d[1]
  Q <- d[2]
  R <- d[3]

  tab <- expand.grid(
    P = seq_len(P),
    Q = seq_len(Q),
    R = seq_len(R)
  )

  tab$Energy <- mapply(
    function(p, q, r) sum(G[p, q, r, ]^2),
    tab$P, tab$Q, tab$R
  )

  tab$RelativeEnergy <- tab$Energy / sum(tab$Energy)
  tab %>% arrange(dplyr::desc(RelativeEnergy))
}

