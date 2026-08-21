test_that("only the two public solvers are exported", {
  expect_setequal(
    getNamespaceExports("fastllqr"),
    c("llqr_seq_ppro", "tvcqr_seq_ppro")
  )
})

test_that("LLQR validates public inputs", {
  x <- seq(0, 1, length.out = 10)
  y <- seq_along(x)
  expect_error(llqr_seq_ppro(x[-1], y), "same length")
  expect_error(llqr_seq_ppro(x, y, tau = 0), "strictly between")
  expect_error(llqr_seq_ppro(x, y, h = 0), "positive finite")
  expect_error(llqr_seq_ppro(x, y, threshold_factor = 0), "positive finite")
  expect_error(llqr_seq_ppro(replace(x, 1, NA), y), "finite")
})

test_that("TVCQR validates public inputs", {
  x <- matrix(seq_len(20), nrow = 10, ncol = 2)
  y <- seq_len(10)
  expect_error(tvcqr_seq_ppro(x, y[-1]), "nrow")
  expect_error(tvcqr_seq_ppro(x, y, tau = 1), "strictly between")
  expect_error(tvcqr_seq_ppro(x, y, h = -1), "positive finite")
  expect_error(tvcqr_seq_ppro(x, y, threshold_factor = Inf), "positive finite")
  expect_error(tvcqr_seq_ppro(replace(x, 1, NaN), y), "finite numeric matrix")
})

test_that("rank-deficient local problems fail explicitly", {
  x <- seq(0, 1, length.out = 20)
  y <- x + stats::rnorm(20)
  expect_error(
    llqr_seq_ppro(x, y, z = 10, h = 0.01, kernel = "epanechnikov"),
    "failed|rank deficient|singular"
  )

  x_tv <- cbind(seq_len(20), seq_len(20))
  expect_error(
    tvcqr_seq_ppro(x_tv, y, h = 0.2),
    "failed|rank deficient|singular|basis"
  )
})
