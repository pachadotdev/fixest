test_that("auto-generated test #73", {
  do_refactor <- FALSE
  do_sum_y <- FALSE
  nthreads <- 4
  obs2keep <- 0L
  only_slope <- FALSE
  r_x_sizes <- 0L
  rm_0 <- FALSE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- structure(list(species = c("setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa"), fe2 = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d")), row.names = c(NA, 49L), class = "data.frame")
  y <- c(5.1, 4.9, 4.7, 4.6, 5, 5.4, 4.6, 5, 4.4, 4.9, 5.4, 4.8, 4.8, 4.3, 5.8, 5.7, 5.4, 5.1, 5.7, 5.1, 5.4, 5.1, 4.6, 5.1, 4.8, 5, 5, 5.2, 5.2, 4.7, 4.8, 5.4, 5.2, 5.5, 4.9, 5, 5.5, 4.9, 4.4, 5.1, 5, 4.5, 4.4, 5, 5.1, 4.8, 5.1, 4.6, 5.3)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), c(9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L)), items = list(1, c(5, 14, 3, 9, 6, 13, 8, 15, 1, 4, 11, 10, 2, 7, 12)), table = list(49L, c(3L, 3L, 4L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 3L, 3L, 4L, 3L, 3L)), sum_y = list(0L, 0L))

  expect_identical(result, expected_result)
})
