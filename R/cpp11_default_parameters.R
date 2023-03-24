# wrappers to define default arguments

cpp_add_commas <- function(x, r = 1L, whole = TRUE) {
    cpp_add_commas_(x, r, whole)
}

cpp_sparse_products <- function(X, w, y, correct_0w = FALSE, nthreads = 1L) {
    cpp_sparse_products_(X, w, y, correct_0w, as.integer(nthreads))
}

cpp_demean <- function(y, X_raw, r_weights, iterMax, diffMax, r_nb_id_Q, fe_id_list, table_id_I, slope_flag_Q, slope_vars_list, r_init, nthreads = 1L, save_fixef = FALSE) {
    cpp_demean_(y, X_raw, r_weights, iterMax, diffMax, r_nb_id_Q, fe_id_list, table_id_I, slope_flag_Q, slope_vars_list, r_init, as.integer(nthreads), save_fixef)
}

cpp_ssq <- function(x, w = 0) {
    cpp_ssq_(x, w)
}



cpppar_crossprod <- function(X, w, nthreads = 1L) {
    cpppar_crossprod_(X, w, as.integer(nthreads))
}

cpppar_xbeta <- function(X, beta, nthreads = 1L) {
    cpppar_xbeta_(X, beta, as.integer(nthreads))
}

cpppar_xwy <- function(X, y, w, nthreads = 1L) {
    cpppar_xwy_(X, y, w, as.integer(nthreads))
}

cpp_tapply_vsum <- function(Q, x, dum) {
    cpp_tapply_vsum_(as.integer(Q), x, as.integer(dum))
}

cpp_partialDerivative_other <- function(iterMax, Q, N, epsDeriv, ll_d2,
                                         dx_dother, init, dumMat, nbCluster) {
    cpp_partialDerivative_other_(as.integer(iterMax), as.integer(Q),
                                 as.integer(N), epsDeriv, ll_d2, dx_dother,
                                 init, dumMat, nbCluster)
}

cpp_conv_seq_poi_2 <- function(n_i, n_j, n_cells, index_i, index_j,
                               dum_vector, sum_y_vector,
                               iterMax, diffMax, exp_mu_in, order) {
    cpp_conv_seq_poi_2_(as.integer(n_i), as.integer(n_j), as.integer(n_cells),
                        index_i, index_j, dum_vector, sum_y_vector,
                        as.integer(iterMax), diffMax, exp_mu_in, order)
}

cpp_conv_seq_gau_2 <- function() {
    cpp_conv_seq_gau_2_(as.integer(n_i), as.integer(n_j), as.integer(n_cells),
                   r_mat_row, r_mat_col, r_mat_value_Ab, r_mat_value_Ba,
                   dum_vector, lhs, invTableCluster_vector,
                   as.integer(iterMax), diffMax, mu_in)
}

cpp_conv_acc_gnl <- function(family, iterMax, diffMax, diffMax_NR, theta,
                             nb_cluster_all, lhs, mu_init, dum_vector,
                             tableCluster_vector, sum_y_vector, cumtable_vector,
                             obsCluster_vector, nthreads) {
    cpp_conv_acc_gnl_(as.integer(family), as.integer(iterMax), diffMax,
                      diffMax_NR, theta, nb_cluster_all, lhs, mu_init,
                      dum_vector, tableCluster_vector, sum_y_vector,
                      cumtable_vector, obsCluster_vector, as.integer(nthreads))
}

cpp_conv_seq_gnl <- function(family, iterMax, diffMax, diffMax_NR, theta,
                             nb_cluster_all, lhs, mu_init, dum_vector,
                             tableCluster_vector, sum_y_vector, cumtable_vector,
                             obsCluster_vector, nthreads) {
    cpp_conv_seq_gnl_(as.integer(family), as.integer(iterMax), diffMax,
                      diffMax_NR, theta, nb_cluster_all, lhs, mu_init,
                      dum_vector, tableCluster_vector, sum_y_vector,
                      cumtable_vector, obsCluster_vector, as.integer(nthreads))
}

cpp_derivconv_seq_2 <- function(iterMax, diffMax, n_vars, nb_cluster_all,
                                n_cells, index_i, index_j, order, ll_d2,
                                jacob_vector, deriv_init_vector, dum_vector) {
    cpp_derivconv_seq_2_(as.integer(iterMax), diffMax, as.integer(n_vars),
                         nb_cluster_all, as.integer(n_cells), index_i, index_j,
                         order, ll_d2, jacob_vector, deriv_init_vector,
                         dum_vector)
}

cpp_derivconv_seq_gnl <- function(iterMax, diffMax, n_vars, nb_cluster_all,
                                   ll_d2, jacob_vector, deriv_init_vector,
                                   dum_vector) {
    cpp_derivconv_seq_gnl_(as.integer(iterMax), diffMax, as.integer(n_vars),
                           nb_cluster_all, ll_d2, jacob_vector,
                           deriv_init_vector, dum_vector)
}

cpp_derivconv_acc_2 <- function(iterMax, diffMax, n_vars, nb_cluster_all,
                                n_cells, index_i, index_j, ll_d2, order,
                                jacob_vector, deriv_init_vector, dum_vector) {
    cpp_derivconv_acc_2_(as.integer(iterMax), diffMax, as.integer(n_vars),
                         nb_cluster_all, as.integer(n_cells), index_i, index_j,
                         ll_d2, order, jacob_vector, deriv_init_vector,
                         dum_vector)
}

cpp_derivconv_acc_gnl <- function(iterMax, diffMax, n_vars, nb_cluster_all,
                                  ll_d2, jacob_vector, deriv_init_vector,
                                  dum_vector) {
    cpp_derivconv_acc_gnl_(as.integer(iterMax), diffMax, as.integer(n_vars),
                           nb_cluster_all, ll_d2, jacob_vector,
                           deriv_init_vector, dum_vector)
}

cpp_fixed_cost_gaussian <- function(n_i, n_cells, index_i, index_j, order,
                                    invTableCluster_vector, dum_vector) {
    cpp_fixed_cost_gaussian_(as.integer(n_i), as.integer(n_cells), index_i,
                            index_j, order, invTableCluster_vector, dum_vector)
}

cpppar_exp <- function(x, nthreads = 1L) {
    cpppar_exp_(x, as.integer(nthreads))
}

cpppar_log <- function(x, nthreads = 1L) {
    cpppar_log_(x, as.integer(nthreads))
}

cpppar_lgamma <- function(x, nthreads = 1L) {
    cpppar_lgamma_(x, as.integer(nthreads))
}

cpppar_digamma <- function(x, nthreads = 1L) {
    cpppar_digamma_(x, as.integer(nthreads))
}

cpppar_trigamma <- function(x, nthreads = 1L) {
    cpppar_trigamma_(x, as.integer(nthreads))
}

cpppar_log_a_exp <- function(nthreads = 1L, a, mu, exp_mu) {
    cpppar_log_a_exp_(as.integer(nthreads), a, mu, exp_mu)
}
