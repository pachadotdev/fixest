test_that("auto-generated test #36", {
  do_refactor <- FALSE
  do_sum_y <- FALSE
  nthreads <- 4
  obs2keep <- 15:126
  only_slope <- FALSE
  r_x_sizes <- 5L
  rm_0 <- FALSE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- list(c(5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L))
  y <- list(0)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L)), items = list(c(5, 1, 2, 3, 4)), table = list(c(24L, 23L, 21L, 21L, 23L)), sum_y = list(0L))

  expect_identical(result, expected_result)
})
