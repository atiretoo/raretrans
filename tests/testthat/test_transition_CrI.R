context("transition_CrI")

# Shared test fixtures
T_mat <- matrix(c(0.5, 0.3, 0.0,
                  0.2, 0.4, 0.1,
                  0.0, 0.1, 0.7), nrow = 3, ncol = 3)
F_mat <- matrix(c(0.0, 0.0, 1.5,
                  0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0), nrow = 3, ncol = 3)
TF <- list(T = T_mat, F = F_mat)
N  <- c(10, 5, 8)

test_that("transition_CrI returns a data frame with correct structure", {
  result <- transition_CrI(TF, N)
  expect_true(is.data.frame(result))
  expect_named(result, c("from_stage", "to_stage", "mean", "lower", "upper"))
})

test_that("transition_CrI returns correct number of rows", {
  # 3 stages x (3 stages + dead) = 12 rows
  result <- transition_CrI(TF, N)
  expect_equal(nrow(result), 12)
})

test_that("transition_CrI probabilities are between 0 and 1", {
  result <- transition_CrI(TF, N)
  expect_true(all(result$mean  >= 0 & result$mean  <= 1))
  expect_true(all(result$lower >= 0 & result$lower <= 1))
  expect_true(all(result$upper >= 0 & result$upper <= 1))
})

test_that("transition_CrI lower <= mean <= upper", {
  result <- transition_CrI(TF, N)
  expect_true(all(result$lower <= result$mean))
  expect_true(all(result$mean  <= result$upper))
})

test_that("transition_CrI stage_names are used correctly", {
  result <- transition_CrI(TF, N, stage_names = c("plantula", "juvenile", "adult"))
  expect_true(all(c("plantula", "juvenile", "adult") %in% result$from_stage))
  expect_true("dead" %in% result$to_stage)
})

test_that("transition_CrI uses colnames when stage_names is NULL", {
  colnames(T_mat) <- c("s1", "s2", "s3")
  TF2 <- list(T = T_mat, F = F_mat)
  result <- transition_CrI(TF2, N)
  expect_true(all(c("s1", "s2", "s3") %in% result$from_stage))
})

test_that("transition_CrI ci parameter changes interval width", {
  result_95 <- transition_CrI(TF, N, ci = 0.95)
  result_50 <- transition_CrI(TF, N, ci = 0.50)
  width_95 <- result_95$upper - result_95$lower
  width_50 <- result_50$upper - result_50$lower
  expect_true(all(width_95 >= width_50))
})

test_that("transition_CrI throws error on bad ci", {
  expect_error(transition_CrI(TF, N, ci = 0))
  expect_error(transition_CrI(TF, N, ci = 1))
  expect_error(transition_CrI(TF, N, ci = -0.5))
  expect_error(transition_CrI(TF, N, ci = c(0.5, 0.9)))
})

test_that("transition_CrI throws error when stage_names wrong length", {
  expect_error(transition_CrI(TF, N, stage_names = c("a", "b")))
})

test_that("transition_CrI throws error on bad TF", {
  expect_error(transition_CrI(N, N))
  expect_error(transition_CrI(list(T = T_mat), N))
})

# --- size boundary tests ---

test_that("transition_CrI works for a 2x2 matrix", {
  T2 <- matrix(c(0.8, 0.1, 0.1, 0.7), nrow = 2)
  F2 <- matrix(c(0.0, 1.2, 0.0, 0.0), nrow = 2)
  TF2 <- list(T = T2, F = F2)
  N2  <- c(20, 15)
  result <- transition_CrI(TF2, N2)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 6)          # 2 stages x (2 + dead) = 6
  expect_true(all(result$lower <= result$mean))
  expect_true(all(result$mean  <= result$upper))
})

test_that("transition_CrI works for a 5x5 matrix", {
  set.seed(42)
  T5 <- matrix(c(
    0.50, 0.05, 0.00, 0.00, 0.00,
    0.20, 0.60, 0.05, 0.00, 0.00,
    0.00, 0.15, 0.65, 0.05, 0.00,
    0.00, 0.00, 0.10, 0.70, 0.05,
    0.00, 0.00, 0.00, 0.10, 0.80), nrow = 5, byrow = TRUE)
  F5 <- matrix(0, nrow = 5, ncol = 5)
  F5[1, 5] <- 2.0
  TF5 <- list(T = T5, F = F5)
  N5  <- c(50, 40, 30, 20, 10)
  result <- transition_CrI(TF5, N5)
  expect_equal(nrow(result), 30)         # 5 stages x (5 + dead) = 30
  expect_true(all(result$lower <= result$mean))
  expect_true(all(result$mean  <= result$upper))
})
