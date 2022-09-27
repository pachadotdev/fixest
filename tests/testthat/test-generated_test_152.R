test_that("auto-generated test #152", {
  do_refactor <- FALSE
  do_sum_y <- FALSE
  nthreads <- 4L
  obs2keep <- 0L
  only_slope <- FALSE
  r_x_sizes <- 0L
  rm_0 <- FALSE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- structure(list(fe1 = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L), fe2 = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L)), class = "data.frame", row.names = c(NA, -20L))
  y <- c(1.03868639960248, 0.0511622852760517, 1.46313562373734, 2.07661883117431, 0.357534660072599, -1.03634206966998, 0.157202327432149, -0.985674286487391, -1.29036652661973, 2.4513795610463, 0.527886904700958, -1.34189750399962, -1.58096732669313, -0.938933220484456, 0.427535629488135, 0.740400921292133, 1.24438381360193, -1.32143423677645, 1.6739874002091, -1.51688470378423)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L), c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L)), items = list(c(1, 2, 3), c(1, 2, 3, 4, 5)), table = list(c(8L, 6L, 6L), c(4L, 4L, 4L, 4L, 4L)), sum_y = list(0L, 0L))

  expect_identical(result, expected_result)
})
