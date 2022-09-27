test_that("auto-generated test #38", {
  do_refactor <- FALSE
  do_sum_y <- FALSE
  nthreads <- 4L
  obs2keep <- 15:126
  only_slope <- FALSE
  r_x_sizes <- 25L
  rm_0 <- FALSE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- list(c(12L, 1L, 2L, 7L, 11L, 13L, 14L, 10L, 3L, 11L, 13L, 9L, 2L, 3L, 11L, 8L, 1L, 2L, 7L, 4L, 8L, 9L, 10L, 3L, 11L, 12L, 14L, 15L, 5L, 16L, 2L, 7L, 4L, 8L, 16L, 17L, 18L, 19L, 20L, 21L, 10L, 18L, 22L, 8L, 14L, 10L, 23L, 4L, 5L, 21L, 10L, 3L, 15L, 20L, 16L, 24L, 7L, 22L, 5L, 9L, 17L, 3L, 22L, 20L, 16L, 2L, 25L, 22L, 20L, 9L, 6L, 3L, 22L, 12L, 16L, 4L, 13L, 9L, 17L, 18L, 4L, 5L, 2L, 7L, 22L, 8L, 14L, 17L, 23L, 22L, 8L, 9L, 10L, 23L, 22L, 13L, 9L, 10L, 18L, 4L, 12L, 16L, 17L, 25L, 4L, 13L, 1L, 17L, 3L, 4L, 8L, 14L))
  y <- list(0)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 5L, 6L, 10L, 3L, 9L, 5L, 11L, 2L, 3L, 4L, 12L, 11L, 10L, 8L, 9L, 5L, 1L, 7L, 13L, 14L, 15L, 3L, 4L, 12L, 11L, 15L, 16L, 17L, 18L, 19L, 20L, 8L, 17L, 21L, 11L, 7L, 8L, 22L, 12L, 14L, 20L, 8L, 9L, 13L, 19L, 15L, 23L, 4L, 21L, 14L, 10L, 16L, 9L, 21L, 19L, 15L, 3L, 24L, 21L, 19L, 10L, 25L, 9L, 21L, 1L, 15L, 12L, 6L, 10L, 16L, 17L, 12L, 14L, 3L, 4L, 21L, 11L, 7L, 16L, 22L, 21L, 11L, 10L, 8L, 22L, 21L, 6L, 10L, 8L, 17L, 12L, 1L, 15L, 16L, 24L, 12L, 6L, 2L, 16L, 9L, 12L, 11L, 7L)), items = list(c(12, 1, 2, 7, 11, 13, 14, 10, 3, 9, 8, 4, 15, 5, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 6)), table = list(c(4L, 3L, 6L, 5L, 4L, 5L, 5L, 7L, 7L, 7L, 7L, 8L, 2L, 4L, 6L, 6L, 4L, 1L, 4L, 2L, 8L, 3L, 1L, 2L, 1L)), sum_y = list(0L))

  expect_identical(result, expected_result)
})
