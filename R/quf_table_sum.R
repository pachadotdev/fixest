# If, for example, you pass in a double to `nthreads` (which would be so easy
# to type `nthreads = 0` instead of `nthreads = 0L`), it'll cause an error.
# So we automatically converts numeric inputs into integers before calling
# `cpppar_quf_table_sum()`, so it doesn't matter if you pass in an integer or a
# double:
quf_table_sum <- function(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope,
                          nthreads, do_refactor, r_x_sizes, obs2keep) {
  cpppar_quf_table_sum(
    x = x, y = y, do_sum_y = do_sum_y, rm_0 = rm_0,
    rm_1 = rm_1, rm_single = rm_single,
    only_slope = only_slope,
    nthreads = as.integer(nthreads),
    do_refactor = do_refactor,
    r_x_sizes = r_x_sizes, obs2keep = as.integer(obs2keep)
  )
}
