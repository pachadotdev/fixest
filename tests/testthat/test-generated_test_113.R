test_that("auto-generated test #113", {
  do_refactor <- TRUE
  do_sum_y <- FALSE
  nthreads <- 4L
  obs2keep <- 0L
  only_slope <- FALSE
  r_x_sizes <- c(1L, 15L)
  rm_0 <- TRUE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- list(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), c(14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L))
  y <- c(6.4, 6.9, 5.5, 6.5, 5.7, 6.3, 4.9, 6.6, 5.2, 5, 5.9, 6, 6.1, 5.6, 6.7, 5.6, 5.8, 6.2, 5.6, 5.9, 6.1, 6.3, 6.1, 6.4, 6.6, 6.8, 6.7, 6, 5.7, 5.5, 5.5, 5.8, 6, 5.4, 6, 6.7, 6.3, 5.6, 5.5, 5.5, 6.1, 5.8, 5, 5.6, 5.7, 5.7, 6.2, 5.1, 5.7)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), c(14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L)), items = list(1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)), table = list(49L, c(3L, 3L, 3L, 4L, 3L, 3L, 4L, 3L, 3L, 3L, 3L, 4L, 3L, 4L, 3L)), sum_y = list(0L, 0L))

  expect_identical(result, expected_result)
})
