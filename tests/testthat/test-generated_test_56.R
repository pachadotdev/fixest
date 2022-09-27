test_that("auto-generated test #56", {
  do_refactor <- FALSE
  do_sum_y <- FALSE
  nthreads <- 4
  obs2keep <- 0
  only_slope <- FALSE
  r_x_sizes <- 0L
  rm_0 <- TRUE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- structure(list(fe_bis = c("b", "k", "y", "o", "k", "y", "i", "u", "g", "w", "s", "c", "f", "z", "o", "l", "p", "z", "e", "i", "v", "j", "u", "z", "z", "s", "q", "z", "v", "x", "m", "o", "w", "k", "v", "n", "r", "p", "y", "j", "j", "t", "x", "h", "t", "t", "m", "k", "h")), row.names = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 17L, 18L, 19L, 20L, 21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L, 29L, 30L, 31L, 32L, 33L, 34L, 35L, 36L, 37L, 38L, 39L, 40L, 41L, 42L, 43L, 44L, 45L, 46L, 47L, 48L, 49L, 50L), class = "data.frame")
  y <- c(5.1, 4.9, 4.7, 4.6, 5, 5.4, 4.6, 5, 4.4, 4.9, 5.4, 4.8, 4.8, 4.3, 5.8, 5.4, 5.1, 5.7, 5.1, 5.4, 5.1, 4.6, 5.1, 4.8, 5, 5, 5.2, 5.2, 4.7, 4.8, 5.4, 5.2, 5.5, 4.9, 5, 5.5, 4.9, 4.4, 5.1, 5, 4.5, 4.4, 5, 5.1, 4.8, 5.1, 4.6, 5.3, 5)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(22L, 19L, 11L, 18L, 19L, 11L, 9L, 20L, 23L, 7L, 17L, 4L, 12L, 10L, 18L, 24L, 8L, 10L, 1L, 9L, 16L, 21L, 20L, 10L, 10L, 17L, 6L, 10L, 16L, 3L, 13L, 18L, 7L, 19L, 16L, 2L, 14L, 8L, 11L, 21L, 21L, 5L, 3L, 15L, 5L, 5L, 13L, 19L, 15L)), items = list(c(19, 36, 30, 12, 42, 27, 10, 17, 7, 14, 3, 13, 31, 37, 44, 21, 11, 4, 2, 8, 22, 1, 9, 16)), table = list(c(1L, 1L, 2L, 1L, 3L, 1L, 2L, 2L, 2L, 5L, 3L, 1L, 2L, 1L, 2L, 3L, 2L, 3L, 4L, 2L, 3L, 1L, 1L, 1L)), sum_y = list(0L))

  expect_identical(result, expected_result)
})
