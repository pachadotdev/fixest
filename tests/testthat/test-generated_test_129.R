test_that("auto-generated test #129", {
  do_refactor <- FALSE
  do_sum_y <- FALSE
  nthreads <- 4
  obs2keep <- 0
  only_slope <- FALSE
  r_x_sizes <- 0L
  rm_0 <- TRUE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- structure(list(species = c("setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa", "setosa"), fe2 = c("g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d")), row.names = 7:49, class = "data.frame")
  y <- c(12.4079547292717, 11.4603019758282, 13.4344739366529, 11.4112129878803, 12.2676507507926, 12.6817527774637, 11.1073046444962, 10.4603310230511, 13.7822289603086, 0, 12.8286045809213, 11.7858067180152, 12.9128286921271, 12.3322651545396, 13.7908192052491, 10.6500741803561, 13.3896382002901, 13.6046516422233, 11.7049377768281, 9.04829361215387, 12.1772373026136, 11.1534418313686, 12.4922032702996, 11.8896367101773, 0, 12.0189604010976, 13.126164353812, 11.8158422795118, 10.7733247260202, 11.0043395013082, 10.024020220854, 10.8974155202363, 10.9409380850995, 11.4534874326665, 11.3664137715686, 9.19089682483287, 10.7582949398644, 13.6535474673476, 12.5224939301207, 0, 11.5945162751585, 11.5091557647992, 11.6658385484487)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), c(13L, 7L, 4L, 11L, 10L, 14L, 6L, 2L, 8L, 12L, 3L, 9L, 1L, 5L, 13L, 7L, 4L, 11L, 10L, 14L, 6L, 2L, 8L, 12L, 3L, 9L, 1L, 5L, 13L, 7L, 4L, 11L, 10L, 14L, 6L, 2L, 8L, 12L, 3L, 9L)), items = list(1, c(14, 8, 12, 3, 15, 7, 2, 9, 13, 5, 4, 11, 1, 6)), table = list(40L, c(2L, 3L, 3L, 3L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L)), sum_y = list(0L, 0L), obs_removed = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE), fe_removed = list(numeric(0), 10))

  expect_identical(result, expected_result)
})
