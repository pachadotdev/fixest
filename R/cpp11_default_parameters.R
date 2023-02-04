# wrappers to define default arguments

cpp_add_commas <- function(x, r = 1L, whole = TRUE) {
    cpp_add_commas_(x, r, whole)
}

cpp_sparse_products <- function(X, w, y, correct_0w = FALSE, nthreads = 1L) {
    cpp_sparse_products_(X, w, y, correct_0w, nthreads)
}

cpp_demean <- function(y, X_raw, r_weights, iterMax, diffMax, r_nb_id_Q, fe_id_list, table_id_I, slope_flag_Q, slope_vars_list, r_init, nthreads, save_fixef = FALSE) {
    cpp_demean_(y, X_raw, r_weights, iterMax, diffMax, r_nb_id_Q, fe_id_list, table_id_I, slope_flag_Q, slope_vars_list, r_init, nthreads, save_fixef)
}

cpp_ssq <- function(x, w = 0) {
    cpp_ssq_(x, w)
}

cpp_partialDerivative_other <- function(iterMax, Q, N, epsDeriv, ll_d2,
                                         dx_dother, init, dumMat, nbCluster) {
    cpp_partialDerivative_other_(iterMax, Q, N, epsDeriv, ll_d2,
                                dx_dother, init, dumMat, nbCluster)
}
