test_that("auto-generated test #80", {
  do_refactor <- FALSE
  do_sum_y <- FALSE
  nthreads <- 4L
  obs2keep <- 0L
  only_slope <- FALSE
  r_x_sizes <- 0L
  rm_0 <- FALSE
  rm_1 <- FALSE
  rm_single <- FALSE
  x <- structure(list(species = c("versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor", "versicolor"), fe2 = c("g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j")), row.names = 52:100, class = "data.frame")
  y <- c(13.6663013557312, 10.9877529466575, 12.4257155118509, 10.352427373005, 9.43412175855268, 11.32702890624, 12.1358625268086, 12.5892298027341, 13.021618767326, 0, 13.5312425187552, 10.3968556670061, 11.6081647625601, 11.9562347974807, 10.7300048943838, 9.50115300503588, 10.8707074085484, 11.1841799043222, 10.3545133893367, 10.6787243344923, 11.7304495034653, 11.1083391910639, 11.8348477615027, 11.3962773739828, 0, 12.7035122324223, 12.2717897755399, 12.5025255953076, 9.88996165933629, 12.1959845897777, 9.50423509662097, 10.8166278574528, 9.97773054885692, 9.29208022119536, 13.5221225187536, 10.8966065891812, 10.8653187806449, 11.1130503963551, 11.6366949746461, 0, 13.1811549557668, 10.1163936117727, 9.7915427464368, 9.837329205279, 10.2468951006308, 13.4093570771456, 11.4576458721328, 10.4073848024104, 10.7988398948477)
  result <- cpppar_quf_table_sum(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads, do_refactor, r_x_sizes, obs2keep)
  expected_result <- list(quf = list(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), c(14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L, 11L, 15L, 6L, 2L, 8L, 9L, 13L, 3L, 10L, 1L, 5L, 14L, 7L, 4L, 12L)), items = list(1, c(14, 8, 12, 3, 15, 7, 2, 9, 10, 13, 5, 4, 11, 1, 6)), table = list(49L, c(3L, 3L, 3L, 4L, 3L, 3L, 4L, 3L, 3L, 3L, 3L, 4L, 3L, 4L, 3L)), sum_y = list(0L, 0L))

  expect_identical_unordered(result, expected_result)
})
