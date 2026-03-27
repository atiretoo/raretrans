context("plot_transition_density")

# Shared test fixtures
T_mat <- matrix(c(0.5, 0.3, 0.0,
                  0.2, 0.4, 0.1,
                  0.0, 0.1, 0.7), nrow = 3, ncol = 3)
F_mat <- matrix(c(0.0, 0.0, 1.5,
                  0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0), nrow = 3, ncol = 3)
TF <- list(T = T_mat, F = F_mat)
N  <- c(10, 5, 8)

test_that("plot_transition_density returns a ggplot object", {
  p <- plot_transition_density(TF, N, stage_names = c("plantula", "juvenile", "adult"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_transition_density runs silently", {
  expect_silent(
    plot_transition_density(TF, N, stage_names = c("plantula", "juvenile", "adult"))
  )
})

test_that("plot_transition_density data has correct columns", {
  p <- plot_transition_density(TF, N)
  expect_true(all(c("from_stage", "to_stage", "x", "density", "ci_density") %in% names(p$data)))
})

test_that("plot_transition_density include_dead = FALSE removes dead row", {
  p_with    <- plot_transition_density(TF, N, include_dead = TRUE)
  p_without <- plot_transition_density(TF, N, include_dead = FALSE)
  expect_gt(nrow(p_with$data), nrow(p_without$data))
  expect_false("dead" %in% as.character(p_without$data$to_stage))
})

test_that("plot_transition_density ci parameter changes shaded region", {
  p_95 <- plot_transition_density(TF, N, ci = 0.95)
  p_50 <- plot_transition_density(TF, N, ci = 0.50)
  # Wider ci means more non-zero ci_density values
  n_shaded_95 <- sum(p_95$data$ci_density > 0, na.rm = TRUE)
  n_shaded_50 <- sum(p_50$data$ci_density > 0, na.rm = TRUE)
  expect_gte(n_shaded_95, n_shaded_50)
})

test_that("plot_transition_density stage_names are used", {
  p <- plot_transition_density(TF, N, stage_names = c("seed", "juvenile", "adult"))
  expect_true(all(c("seed", "juvenile", "adult") %in% as.character(p$data$from_stage)))
})

test_that("plot_transition_density throws error on bad ci", {
  expect_error(plot_transition_density(TF, N, ci = 0))
  expect_error(plot_transition_density(TF, N, ci = 1))
  expect_error(plot_transition_density(TF, N, ci = 1.5))
})

test_that("plot_transition_density throws error when stage_names wrong length", {
  expect_error(plot_transition_density(TF, N, stage_names = c("a", "b")))
})

test_that("plot_transition_density throws error on bad TF", {
  expect_error(plot_transition_density(N, N))
  expect_error(plot_transition_density(list(T = T_mat), N))
})

# --- size boundary tests ---

test_that("plot_transition_density works for a 2x2 matrix", {
  T2  <- matrix(c(0.8, 0.1, 0.1, 0.7), nrow = 2)
  F2  <- matrix(c(0.0, 1.2, 0.0, 0.0), nrow = 2)
  TF2 <- list(T = T2, F = F2)
  N2  <- c(20, 15)
  p <- plot_transition_density(TF2, N2, stage_names = c("juvenile", "adult"))
  expect_s3_class(p, "ggplot")
  # 2 stages x 3 fates x 300 x-points = 1800 rows
  expect_equal(nrow(p$data), 1800)
})

test_that("plot_transition_density works for a 5x5 matrix", {
  T5 <- matrix(c(
    0.50, 0.05, 0.00, 0.00, 0.00,
    0.20, 0.60, 0.05, 0.00, 0.00,
    0.00, 0.15, 0.65, 0.05, 0.00,
    0.00, 0.00, 0.10, 0.70, 0.05,
    0.00, 0.00, 0.00, 0.10, 0.80), nrow = 5, byrow = TRUE)
  F5 <- matrix(0, nrow = 5, ncol = 5); F5[1, 5] <- 2.0
  TF5 <- list(T = T5, F = F5)
  N5  <- c(50, 40, 30, 20, 10)
  p <- plot_transition_density(TF5, N5,
                               stage_names = c("s1", "s2", "s3", "s4", "s5"))
  expect_s3_class(p, "ggplot")
  # 5 stages x 6 fates x 300 x-points = 9000 rows
  expect_equal(nrow(p$data), 9000)
})
