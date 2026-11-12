test_that("Fourier basis has expected dimensions and finite values", {
  Phi <- fourier_basis_matrix(0:7, N = 3, L = 8)
  expect_equal(dim(Phi), c(8L, 3L))
  expect_true(all(is.finite(Phi)))
})

test_that("Fourier coefficients and centering can be reconstructed", {
  set.seed(22)
  X <- array(rnorm(3 * 4 * 3 * 8), c(3, 4, 3, 8))
  cf <- build_coefficient_tensor(X, 0:7, N = 3, L = 8)
  expect_equal(dim(cf$C), c(3L, 4L, 3L, 3L))
  centered <- center_coefficient_tensor(cf$C)
  expect_equal(decenter_coefficient_tensor(centered$C_centered, centered),
               cf$C, tolerance = 1e-9)
})

test_that("FTUCKALS-3 uses three reduced modes and preserves basis mode", {
  set.seed(23)
  C <- array(rnorm(4 * 4 * 3 * 3), c(4, 4, 3, 3))
  fit <- ftuckals3(C, P = 2, Q = 2, R = 1, max_iter = 15)
  expect_equal(dim(fit$G), c(2L, 2L, 1L, 3L))
  expect_equal(dim(fit$C_hat), dim(C))
  expect_true(is.finite(fit$relative_error))
})

test_that("the public interface returns reconstructed functional observations", {
  set.seed(24)
  X <- array(rnorm(3 * 4 * 3 * 8), c(3, 4, 3, 8))
  result <- ftucker3(X, 0:7, ranks = c(2, 2, 1), N = 3,
                     L = 8, max_iter = 10)
  expect_s3_class(result, "FTucker3Fit")
  expect_equal(dim(result$reconstructed), dim(X))
  expect_equal(dim(predict(result)), dim(X))
})
