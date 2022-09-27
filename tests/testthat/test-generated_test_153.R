test_that("auto-generated test #153", {
  do_refactor <- FALSE
  do_sum_y <- FALSE
  nthreads <- 4
  obs2keep <- 0
  only_slope <- FALSE
  r_x_sizes <- 0L
  rm_0 <- FALSE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- structure(list(grp = c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L), tm = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L)), class = "data.frame", row.names = c(NA, -20L))
  y <- c(1.26295428488079, -0.326233360705649, 1.3297992629225, 1.2724293214294, 0.414641434456408, -1.53995004190371, -0.928567034713538, -0.29472044679056, -0.00576717274753696, 2.40465338885795, 0.76359346114046, -0.799009248989368, -1.14765700923635, -0.289461573688223, -0.299215117897316, -0.411510832795067, 0.252223448156132, -0.891921127284569, 0.435683299355719, -1.23753842192996)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L), c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L)), items = list(c(1, 2, 3, 4), c(1, 2, 3, 4, 5)), table = list(c(5L, 5L, 5L, 5L), c(4L, 4L, 4L, 4L, 4L)), sum_y = list(0L, 0L))

  expect_identical(result, expected_result)
})
