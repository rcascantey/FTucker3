test_that("FTUCKALS-3 reconstructs an exact low-rank tensor", {
  
  set.seed(20261112)
  
  I <- 5
  J <- 4
  K <- 3
  N <- 3
  
  P <- 2
  Q <- 2
  R <- 1
  
  # Structural loading matrices
  A <- qr.Q(qr(matrix(rnorm(I * P), I, P)))
  B <- qr.Q(qr(matrix(rnorm(J * Q), J, Q)))
  C <- qr.Q(qr(matrix(rnorm(K * R), K, R)))
  
  # Functional coefficient core
  G <- array(
    rnorm(P * Q * R * N),
    dim = c(P, Q, R, N)
  )
  
  # Exact coefficient tensor
  Ctrue <- array(0, dim = c(I, J, K, N))
  
  for (n in seq_len(N)) {
    for (k in seq_len(K)) {
      Ctrue[, , k, n] <-
        C[k, 1] * (
          A %*% G[, , 1, n] %*% t(B)
        )
    }
  }
  
  # Estimate FTucker-3
  fit <- ftuckals3(
    Ccenter = Ctrue,
    P = P,
    Q = Q,
    R = R,
    max_iter = 200,
    tol = 1e-10,
    n_random_starts = 0
  )
  
  # Relative reconstruction error
  error <- sqrt(
    sum((Ctrue - fit$C_hat)^2) /
      sum(Ctrue^2)
  )
  
  # Automatic checks
  expect_equal(
    dim(fit$G),
    c(P, Q, R, N)
  )
  
  expect_lt(error, 1e-7)
  
  expect_equal(
    fit$R2_coeff,
    1,
    tolerance = 1e-7
  )
  
  expect_true(fit$converged)
  
})
