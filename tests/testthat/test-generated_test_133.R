test_that("auto-generated test #133", {
  do_refactor <- FALSE
  do_sum_y <- FALSE
  nthreads <- 4
  obs2keep <- 0
  only_slope <- FALSE
  r_x_sizes <- 0L
  rm_0 <- TRUE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- structure(list(species = c("virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica", "virginica"), fe2 = c("k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o")), row.names = 101:150, class = "data.frame")
  y <- c(12.7244594064128, 11.6105978350922, 11.1857280198281, 10.7003699045169, 10.6378016703142, 0, 12.1899200776134, 13.4586871159153, 10.8286264275946, 11.4491655768521, 10.5726194019132, 11.0994808732797, 11.9718594661169, 12.6089398209981, 11.9641686026836, 12.0559800904812, 12.7309536630209, 13.0471368477271, 11.4065980409274, 10.3166833429914, 0, 11.5388584191035, 13.1106315882366, 10.919359025277, 10.6057704190729, 12.1375795245805, 10.7304140128498, 12.1388056114383, 9.67601016550743, 9.75756991965909, 12.089804173283, 12.2309631771735, 12.2710677089481, 9.38375441778656, 12.5125791035104, 0, 12.7320683259354, 12.3364102561772, 13.6100735137793, 10.0961901531867, 10.9668961518694, 11.9597239825503, 10.5430183645858, 11.6855504408545, 12.3962431686397, 10.8463269386689, 11.907105983302, 12.0499092350097, 10.8932706415673, 10.5026202827238)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), c(10L, 14L, 6L, 2L, 8L, 12L, 3L, 9L, 1L, 5L, 13L, 7L, 4L, 11L, 10L, 14L, 6L, 2L, 8L, 12L, 3L, 9L, 1L, 5L, 13L, 7L, 4L, 11L, 10L, 14L, 6L, 2L, 8L, 12L, 3L, 9L, 1L, 5L, 13L, 7L, 4L, 11L, 10L, 14L, 6L, 2L, 8L)), items = list(1, c(10, 4, 8, 14, 11, 3, 13, 5, 9, 1, 15, 7, 12, 2)), table = list(47L, c(3L, 4L, 3L, 3L, 3L, 4L, 3L, 4L, 3L, 4L, 3L, 3L, 3L, 4L)), sum_y = list(0L, 0L), obs_removed = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), fe_removed = list(numeric(0), 6))

  expect_identical(result, expected_result)
})
