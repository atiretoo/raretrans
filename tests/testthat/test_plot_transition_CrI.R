context("plot_transition_CrI")

# Shared test fixtures
T_mat <- matrix(c(0.5, 0.3, 0.0,
                  0.2, 0.4, 0.1,
                  0.0, 0.1, 0.7), nrow = 3, ncol = 3)
F_mat <- matrix(c(0.0, 0.0, 1.5,
                  0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0), nrow = 3, ncol = 3)
TF  <- list(T = T_mat, F = F_mat)
N   <- c(10, 5, 8)
cri <- transition_CrI(TF, N, stage_names = c("plantula", "juvenile", "adult"))

test_that("plot_transition_CrI returns a ggplot object", {
  p <- plot_transition_CrI(cri)
  expect_s3_class(p, "ggplot")
})

test_that("plot_transition_CrI runs silently", {
  expect_silent(plot_transition_CrI(cri))
})

test_that("plot_transition_CrI include_dead = FALSE removes dead row", {
  p_with    <- plot_transition_CrI(cri, include_dead = TRUE)
  p_without <- plot_transition_CrI(cri, include_dead = FALSE)
  # The data underlying the plot should differ in number of rows
  expect_gt(nrow(p_with$data),    nrow(p_without$data))
  expect_false("dead" %in% p_without$data$to_stage)
})

test_that("plot_transition_CrI title can be changed", {
  p <- plot_transition_CrI(cri, title = "My custom title")
  expect_equal(p$labels$title, "My custom title")
})

test_that("plot_transition_CrI throws error on bad cri input", {
  expect_error(plot_transition_CrI(data.frame(a = 1, b = 2)))
  expect_error(plot_transition_CrI("not a data frame"))
  expect_error(plot_transition_CrI(list(from_stage = 1)))
})

# --- size boundary tests ---

test_that("plot_transition_CrI works for a 2x2 matrix", {
  T2  <- matrix(c(0.8, 0.1, 0.1, 0.7), nrow = 2)
  F2  <- matrix(c(0.0, 1.2, 0.0, 0.0), nrow = 2)
  TF2 <- list(T = T2, F = F2)
  N2  <- c(20, 15)
  cri2 <- transition_CrI(TF2, N2, stage_names = c("juvenile", "adult"))
  p <- plot_transition_CrI(cri2)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 6)   # 2 stages x 3 fates (incl. dead)
})

test_that("plot_transition_CrI works for a 5x5 matrix", {
  T5 <- matrix(c(
    0.50, 0.05, 0.00, 0.00, 0.00,
    0.20, 0.60, 0.05, 0.00, 0.00,
    0.00, 0.15, 0.65, 0.05, 0.00,
    0.00, 0.00, 0.10, 0.70, 0.05,
    0.00, 0.00, 0.00, 0.10, 0.80), nrow = 5, byrow = TRUE)
  F5 <- matrix(0, nrow = 5, ncol = 5); F5[1, 5] <- 2.0
  TF5 <- list(T = T5, F = F5)
  N5  <- c(50, 40, 30, 20, 10)
  cri5 <- transition_CrI(TF5, N5,
                         stage_names = c("s1", "s2", "s3", "s4", "s5"))
  p <- plot_transition_CrI(cri5)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 30)  # 5 stages x 6 fates (incl. dead)
})
