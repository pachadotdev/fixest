# If cpp11.R is changed, this file needs to be updated.

cpp_fixed_cost_gaussian <- function(n_i, n_cells, index_i, index_j, order, invTableCluster_vector, dum_vector) {
  cpp_fixed_cost_gaussian_(parse_int(n_i), parse_int(n_cells), parse_int(index_i), parse_int(index_j), order, invTableCluster_vector, dum_vector)
}

update_mu_single_cluster <- function(family, nb_cluster, theta, diffMax_NR, mu_in, lhs, sum_y, dum, obsCluster, table, cumtable, nthreads = 1L) {
  update_mu_single_cluster_(parse_int(family), parse_int(nb_cluster), theta, diffMax_NR, mu_in, lhs, sum_y, parse_int(dum), parse_int_mat(obsCluster), table, cumtable, parse_int(nthreads))
}

get_n_cells <- function(index_i, index_j) {
  get_n_cells_(parse_int(index_i), parse_int(index_j))
}

cpp_derivconv_seq_gnl <- function(iterMax, diffMax, n_vars, nb_cluster_all, ll_d2, jacob_vector, deriv_init_vector, dum_vector) {
  cpp_derivconv_seq_gnl_(parse_int(iterMax), diffMax, n_vars, nb_cluster_all, ll_d2, jacob_vector, deriv_init_vector, dum_vector)
}

cpp_derivconv_acc_gnl <- function(iterMax, diffMax, n_vars, nb_cluster_all, ll_d2, jacob_vector, deriv_init_vector, dum_vector) {
  cpp_derivconv_acc_gnl_(parse_int(iterMax), diffMax, n_vars, nb_cluster_all, ll_d2, jacob_vector, deriv_init_vector, dum_vector)
}

cpp_derivconv_acc_2 <- function(iterMax, diffMax, n_vars, nb_cluster_all, n_cells, index_i, index_j, ll_d2, order, jacob_vector, deriv_init_vector, dum_vector) {
  cpp_derivconv_acc_2_(parse_int(iterMax), diffMax, n_vars, nb_cluster_all, parse_int(n_cells), parse_int(index_i), parse_int(index_j), ll_d2, order, jacob_vector, deriv_init_vector, dum_vector)
}

cpp_derivconv_seq_2 <- function(iterMax, diffMax, n_vars, nb_cluster_all, n_cells, index_i, index_j, order, ll_d2, jacob_vector, deriv_init_vector, dum_vector) {
  cpp_derivconv_seq_2_(parse_int(iterMax), diffMax, n_vars, nb_cluster_all, parse_int(n_cells), parse_int(index_i), parse_int(index_j), order, ll_d2, jacob_vector, deriv_init_vector, dum_vector)
}

update_deriv_single <- function(n_vars, nb_coef, r_ll_d2, r_jacob_vector, r_dum_vector) {
  update_deriv_single_(n_vars, parse_int(nb_coef), r_ll_d2, r_jacob_vector, r_dum_vector)
}

cpp_conv_acc_gnl <- function(family, iterMax, diffMax, diffMax_NR, theta, nb_cluster_all, lhs, mu_init, dum_vector, tableCluster_vector, sum_y_vector, cumtable_vector, obsCluster_vector, nthreads = 1L) {
  cpp_conv_acc_gnl_(parse_int(family), parse_int(iterMax), diffMax, diffMax_NR, theta, nb_cluster_all, lhs, mu_init, dum_vector, tableCluster_vector, sum_y_vector, cumtable_vector, obsCluster_vector, parse_int(nthreads))
}

cpp_conv_seq_gnl <- function(family, iterMax, diffMax, diffMax_NR, theta, nb_cluster_all, lhs, mu_init, dum_vector, tableCluster_vector, sum_y_vector, cumtable_vector, obsCluster_vector, nthreads = 1L) {
  cpp_conv_seq_gnl_(parse_int(family), parse_int(iterMax), diffMax, diffMax_NR, theta, nb_cluster_all, lhs, mu_init, dum_vector, tableCluster_vector, sum_y_vector, cumtable_vector, obsCluster_vector, parse_int(nthreads))
}

cpp_conv_acc_poi_2 <- function(n_i, n_j, n_cells, index_i, index_j, dum_vector, sum_y_vector, iterMax, diffMax, exp_mu_in, order) {
  cpp_conv_acc_poi_2_(parse_int(n_i), parse_int(n_j), parse_int(n_cells), parse_int(index_i), parse_int(index_j), dum_vector, sum_y_vector, parse_int(iterMax), diffMax, exp_mu_in, order)
}

cpp_conv_seq_poi_2 <- function(n_i, n_j, n_cells, index_i, index_j, dum_vector, sum_y_vector, iterMax, diffMax, exp_mu_in, order) {
  cpp_conv_seq_poi_2_(parse_int(n_i), parse_int(n_j), parse_int(n_cells), parse_int(index_i), parse_int(index_j), dum_vector, sum_y_vector, parse_int(iterMax), diffMax, exp_mu_in, order)
}

cpp_conv_acc_gau_2 <- function(n_i, n_j, n_cells, r_mat_row, r_mat_col, r_mat_value_Ab, r_mat_value_Ba, dum_vector, lhs, invTableCluster_vector, iterMax, diffMax, mu_in) {
  cpp_conv_acc_gau_2_(parse_int(n_i), parse_int(n_j), parse_int(n_cells), r_mat_row, r_mat_col, r_mat_value_Ab, r_mat_value_Ba, dum_vector, lhs, invTableCluster_vector, parse_int(iterMax), diffMax, mu_in)
}

cpp_conv_seq_gau_2 <- function(n_i, n_j, n_cells, r_mat_row, r_mat_col, r_mat_value_Ab, r_mat_value_Ba, dum_vector, lhs, invTableCluster_vector, iterMax, diffMax, mu_in) {
  cpp_conv_seq_gau_2_(parse_int(n_i), parse_int(n_j), parse_int(n_cells), r_mat_row, r_mat_col, r_mat_value_Ab, r_mat_value_Ba, dum_vector, lhs, invTableCluster_vector, parse_int(iterMax), diffMax, mu_in)
}

cpp_which_na_inf <- function(x, nthreads = 1L) {
  cpp_which_na_inf_(x, parse_int(nthreads))
}

cpp_demean <- function(y, X_raw, r_weights, iterMax, diffMax, r_nb_id_Q, fe_id_list, table_id_I, slope_flag_Q, slope_vars_list, r_init, nthreads = 1L, algo_extraProj, algo_iter_warmup, algo_iter_projAfterAcc, algo_iter_grandAcc, save_fixef = FALSE) {
  cpp_demean_(y, X_raw, r_weights, parse_int(iterMax), diffMax, r_nb_id_Q, fe_id_list, table_id_I, slope_flag_Q, slope_vars_list, r_init, parse_int(nthreads), algo_extraProj, algo_iter_warmup, algo_iter_projAfterAcc, algo_iter_grandAcc, save_fixef)
}

cpp_dsb <- function(Rstr) {
  cpp_dsb_(Rstr)
}

cpp_dsb_full_string <- function(Rstr) {
  cpp_dsb_full_string_(Rstr)
}

cpp_dsb_if_extract <- function(Rstr) {
  cpp_dsb_if_extract_(Rstr)
}

cpp_paste_conditional <- function(x, id, n) {
  cpp_paste_conditional_(x, parse_int(id), parse_int(n))
}

cpp_cholesky <- function(X, tol = 1e-10, nthreads = 1L) {
  cpp_cholesky_(X, tol, parse_int(nthreads))
}

cpp_sparse_products <- function(X, w, y, correct_0w = FALSE, nthreads = 1L) {
  cpp_sparse_products_(X, w, y, correct_0w, parse_int(nthreads))
}

cpp_crossprod <- function(X, w, nthreads = 1L) {
  cpp_crossprod_(X, w, parse_int(nthreads))
}

cpp_mat_reconstruct <- function(X, id_excl) {
  cpp_mat_reconstruct_(X, id_excl)
}

cpp_iv_products <- function(X, y, Z, u, w, nthreads = 1L) {
  cpp_iv_products_(X, y, Z, u, w, parse_int(nthreads))
}

cpp_iv_product_completion <- function(XtX, Xty, X, y, U, w, nthreads = 1L) {
  cpp_iv_product_completion_(XtX, Xty, X, y, U, w, parse_int(nthreads))
}

cpp_iv_resid <- function(resid_2nd, coef, resid_1st, is_int, nthreads = 1L) {
  cpp_iv_resid_(resid_2nd, coef, resid_1st, is_int, parse_int(nthreads))
}

cpp_add_commas <- function(x, r = 1L, whole = TRUE) {
  cpp_add_commas_(x, r, whole)
}

cpp_find_never_always_treated <- function(cohort, period) {
  cpp_find_never_always_treated_(cohort, period)
}

cpp_get_first_item <- function(x, n_items) {
  cpp_get_first_item_(x, n_items)
}

cpp_combine_clusters <- function(cluster_list, index) {
  cpp_combine_clusters_(cluster_list, index)
}

cpp_cut <- function(x_sorted, cut_points, is_included) {
  cpp_cut_(x_sorted, cut_points, is_included)
}

cpp_is_int <- function(x) {
  cpp_is_int_(x)
}

cpp_hash_string <- function(x) {
  cpp_hash_string_(x)
}

cpp_isConstant <- function(x) {
  cpp_isConstant_(x)
}

cpp_any_na_null <- function(x) {
  cpp_any_na_null_(x)
}

cpp_find_duplicates <- function(id, time) {
  cpp_find_duplicates_(parse_int(id), parse_int(time))
}

cpp_check_nested <- function(fe_list, cluster_list, fe_sizes, n) {
  cpp_check_nested_(fe_list, cluster_list, parse_int_mat(fe_sizes), parse_int(n))
}

cpp_pgcd <- function(x) {
  cpp_pgcd_(parse_int(x))
}

cpp_escape_markup <- function(Rstr) {
  cpp_escape_markup_(Rstr)
}

cpp_factor_matrix <- function(fact, is_na_all, who_is_dropped, var, col_names) {
  cpp_factor_matrix_(parse_int(fact), is_na_all, parse_int(who_is_dropped), var, parse_chr(col_names))
}

cpp_diag_XUtX <- function(X, U) {
  cpp_diag_XUtX_(X, U)
}

cpp_partialDerivative_other <- function(iterMax, Q, N, epsDeriv, ll_d2, dx_dother, init, dumMat, nbCluster) {
  cpp_partialDerivative_other_(parse_int(iterMax), parse_int(Q), parse_int(N), epsDeriv, ll_d2, dx_dother, init, parse_int_mat(dumMat), nbCluster)
}

cpp_table <- function(Q, dum) {
  cpp_table_(parse_int(Q), parse_int(dum))
}

cpp_ssr_null <- function(y, w = 0) {
  cpp_ssr_null_(y, w)
}

cpp_ssq <- function(x, w = 0) {
  cpp_ssq_(x, w)
}

cpp_constant_dum <- function(k, x, dum, only_0 = FALSE) {
  cpp_constant_dum_(k, x, parse_int(dum), only_0)
}

cpp_tapply_sum <- function(Q, x, dum) {
  cpp_tapply_sum_(parse_int(Q), x, parse_int(dum))
}

cpp_tapply_vsum <- function(Q, x, dum) {
  cpp_tapply_vsum_(parse_int(Q), x, parse_int(dum))
}

cpp_lag_obs <- function(id, time, nlag) {
  cpp_lag_obs_(parse_int(id), parse_int(time), parse_int(nlag))
}

cpp_get_fe_gnl <- function(Q, N, sumFE, dumMat, cluster_sizes, obsCluster) {
  cpp_get_fe_gnl_(parse_int(Q), parse_int(N), sumFE, parse_int_mat(dumMat), parse_int(cluster_sizes), parse_int_mat(obsCluster))
}

cpp_quf_gnl <- function(x) {
  cpp_quf_gnl_(x)
}

cpp_quf_table_sum <- function(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads = 1L, do_refactor, r_x_sizes, obs2keep) {
  cpp_quf_table_sum_(x, y, do_sum_y, rm_0, rm_1, rm_single, parse_int(only_slope), parse_int(nthreads), do_refactor, parse_int(r_x_sizes), parse_int(obs2keep))
}

cpp_get_nb_threads <- function() {
  cpp_get_nb_threads_()
}

cpp_which_na_inf_vec <- function(x, nthreads = 1L) {
  cpp_which_na_inf_vec_(x, parse_int(nthreads))
}

cpp_which_na_inf_mat <- function(mat, nthreads = 1L) {
  cpp_which_na_inf_mat_(mat, parse_int(nthreads))
}

cpp_which_na_inf_df <- function(df, nthreads = 1L) {
  cpp_which_na_inf_df_(df, parse_int(nthreads))
}

cpp_check_only_0 <- function(x_mat, nthreads = 1L) {
  cpp_check_only_0_(x_mat, parse_int(nthreads))
}

cpp_exp <- function(x, nthreads = 1L) {
  cpp_exp_(x, parse_int(parse_int(nthreads)))
}

cpp_log <- function(x, nthreads = 1L) {
  cpp_log_(x, parse_int(parse_int(nthreads)))
}

cpp_log_a_exp <- function(a, mu, exp_mu, nthreads = 1L) {
  cpp_log_a_exp_(a, mu, exp_mu, parse_int(parse_int(nthreads)))
}

cpp_cond_means <- function(mat_vars, treat, nthreads = 1L) {
  cpp_cond_means_(mat_vars, parse_int(treat), parse_int(nthreads))
}

cpp_lgamma <- function(x, nthreads = 1L) {
  cpp_lgamma_(x, parse_int(parse_int(nthreads)))
}

cpp_digamma <- function(x, nthreads = 1L) {
  cpp_digamma_(x, parse_int(parse_int(nthreads)))
}

cpp_trigamma <- function(x, nthreads = 1L) {
  cpp_trigamma_(x, parse_int(parse_int(nthreads)))
}

cpp_poisson_linkinv <- function(x, nthreads = 1L) {
  cpp_poisson_linkinv_(x, parse_int(nthreads))
}

cpp_poisson_validmu <- function(x, nthreads = 1L) {
  cpp_poisson_validmu_(x, parse_int(parse_int(nthreads)))
}

cpp_logit_linkfun <- function(x, nthreads = 1L) {
  cpp_logit_linkfun_(x, parse_int(parse_int(nthreads)))
}

cpp_logit_linkinv <- function(x, nthreads = 1L) {
  cpp_logit_linkinv_(x, parse_int(parse_int(nthreads)))
}

cpp_logit_mueta <- function(x, nthreads = 1L) {
  cpp_logit_mueta_(x, parse_int(parse_int(nthreads)))
}

cpp_logit_devresids <- function(y, mu, wt, nthreads = 1L) {
  cpp_logit_devresids_(y, mu, wt, parse_int(parse_int(nthreads)))
}

cpp_xwy <- function(X, y, w, nthreads = 1L) {
  cpp_xwy_(X, y, w, parse_int(nthreads))
}

cpp_xbeta <- function(X, beta, nthreads = 1L) {
  cpp_xbeta_(X, beta, parse_int(nthreads))
}

cpp_matprod <- function(x, y, nthreads = 1L) {
  cpp_matprod_(x, y, parse_int(nthreads))
}

cpp_colon_to_star <- function(Rstr) {
  cpp_colon_to_star_(Rstr)
}

cpp_newey_west <- function(S, w, nthreads = 1L) {
  cpp_newey_west_(S, w, parse_int(nthreads))
}

cpp_newey_west_panel <- function(S, w, unit, G, time, T, nthreads = 1L) {
  cpp_newey_west_panel_(S, w, parse_int(unit), parse_int(G), parse_int(time), parse_int(T), parse_int(nthreads))
}

cpp_driscoll_kraay <- function(S, w, time, T, nthreads = 1L) {
  cpp_driscoll_kraay_(S, w, parse_int(time), parse_int(T), parse_int(nthreads))
}

cpp_vcov_conley <- function(S, lon_rad, lat_rad, distance, cutoff, nthreads = 1L) {
  cpp_vcov_conley_(S, lon_rad, lat_rad, parse_int(distance), cutoff, parse_int(nthreads))
}
