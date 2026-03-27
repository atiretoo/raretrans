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
