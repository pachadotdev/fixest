expect_identical_unordered <- function(result, expected_result){
    sorted_result <- purrr::map_depth(result, 2, sort)
    sorted_expected <- purrr::map_depth(expected_result, 2, sort)

    # as the values assigned aren't important; just that they're consistent
    sorted_expected$quf <- purrr::map(sorted_expected$quf, order)
    sorted_result$quf <- purrr::map(sorted_result$quf, order)

    # expected values should be doubles
    sorted_expected$sum_y <- purrr::map(sorted_expected$sum_y, as.double)

    expect_equal(sorted_result, sorted_expected)
}
