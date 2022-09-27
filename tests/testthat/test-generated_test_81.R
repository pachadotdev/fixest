test_that("auto-generated test #81", {
  do_refactor <- FALSE
  do_sum_y <- FALSE
  nthreads <- 4
  obs2keep <- 0L
  only_slope <- FALSE
  r_x_sizes <- 0L
  rm_0 <- FALSE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- structure(list(species = c("virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica"), fe2 = c("k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o")), row.names = 101:150, class = "data.frame")
  y <- c(6.3, 5.8, 7.1, 6.3, 6.5, 7.6, 4.9, 7.3, 6.7, 7.2, 6.5, 6.4, 6.8, 5.7, 5.8, 6.4, 6.5, 7.7, 7.7, 6, 6.9, 5.6, 7.7, 6.3, 6.7, 7.2, 6.2, 6.1, 6.4, 7.2, 7.4, 7.9, 6.4, 6.3, 6.1, 7.7, 6.3, 6.4, 6, 6.9, 6.7, 6.9, 5.8, 6.8, 6.7, 6.7, 6.3, 6.5, 6.2, 5.9)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), c(11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L)), items = list(1, c(10, 4, 8, 14, 11, 3, 13, 5, 6, 9, 1, 15, 7, 12, 2)), table = list(50L, c(3L, 4L, 3L, 3L, 3L, 4L, 3L, 4L, 3L, 3L, 4L, 3L, 3L, 3L, 4L)), sum_y = list(0L, 0L))

  expect_identical(result, expected_result)
})
