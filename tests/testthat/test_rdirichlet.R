context("rdirichlet")

test_that("rdirichlet returns a matrix with correct dimensions", {
  set.seed(42)
  result <- rdirichlet(5, c(1, 2, 3))
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 3)
})

test_that("each row of rdirichlet sums to 1", {
  set.seed(42)
  result <- rdirichlet(10, c(1, 2, 3))
  expect_equal(rowSums(result), rep(1, 10), tolerance = 1e-10)
})

test_that("rdirichlet values are all non-negative", {
  set.seed(42)
  result <- rdirichlet(20, c(0.5, 1, 2, 3))
  expect_true(all(result >= 0))
})

test_that("rdirichlet with n = 1 returns a single-row matrix", {
  set.seed(42)
  result <- rdirichlet(1, c(1, 1, 1))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 3)
  expect_equal(sum(result), 1, tolerance = 1e-10)
})

test_that("rdirichlet is reproducible with set.seed", {
  set.seed(99)
  r1 <- rdirichlet(5, c(1, 2, 3))
  set.seed(99)
  r2 <- rdirichlet(5, c(1, 2, 3))
  expect_equal(r1, r2)
})

test_that("rdirichlet respects concentration: high alpha concentrates near mean", {
  set.seed(42)
  # With very high alpha, samples should be close to the mean (1/3, 1/3, 1/3)
  result <- rdirichlet(1000, c(1000, 1000, 1000))
  col_means <- colMeans(result)
  expect_equal(col_means, c(1/3, 1/3, 1/3), tolerance = 0.01)
})
