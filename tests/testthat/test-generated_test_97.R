test_that("auto-generated test #97", {
  do_refactor <- TRUE
  do_sum_y <- FALSE
  nthreads <- 4L
  obs2keep <- 0L
  only_slope <- FALSE
  r_x_sizes <- c(15L, 1L)
  rm_0 <- FALSE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- list(c(11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L), c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L))
  y <- c(6.3, 5.8, 7.1, 6.3, 6.5, 7.6, 4.9, 7.3, 6.7, 7.2, 6.5, 6.4, 6.8, 5.7, 5.8, 6.4, 6.5, 7.7, 7.7, 6, 6.9, 5.6, 7.7, 6.3, 6.7, 7.2, 6.2, 6.1, 6.4, 7.2, 7.4, 7.9, 6.4, 6.3, 6.1, 7.7, 6.3, 6.4, 6, 6.9, 6.7, 6.9, 5.8, 6.8, 6.7, 6.7, 6.3, 6.5, 6.2, 5.9)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L), c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L)), items = list(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), 1), table = list(c(3L, 4L, 3L, 3L, 3L, 4L, 3L, 4L, 3L, 3L, 4L, 3L, 3L, 3L, 4L), 50L), sum_y = list(0L, 0L))

  expect_identical_unordered(result, expected_result)
})
