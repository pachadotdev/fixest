# This file needs to be manually updated when the C++ code is changed.
# i.e., if cpp11.R is changed, this file needs to be updated.

cpp_fixed_cost_gaussian <- function(n_i, n_cells, index_i, index_j, order, invTableCluster_vector, dum_vector) {
  .Call(`_fixest2_cpp_fixed_cost_gaussian_`, as.integer(n_i), as.integer(n_cells), index_i, index_j, order, invTableCluster_vector, dum_vector)
}

compute_cluster_coef_r <- function(family, nb_coef, theta, diffMax_NR, r_mu, r_lhs, r_sum_y, r_dum, r_obsCluster, r_table, r_cumtable, nthreads = 1L) {
  .Call(`_fixest2_compute_cluster_coef_r_`, as.integer(family), as.integer(nb_coef), theta, diffMax_NR, r_mu, r_lhs, r_sum_y, r_dum, r_obsCluster, r_table, r_cumtable, as.integer(nthreads))
}

update_mu_single_cluster <- function(family, nb_cluster, theta, diffMax_NR, mu_in, lhs, sum_y, dum, obsCluster, table, cumtable, nthreads = 1L) {
  .Call(`_fixest2_update_mu_single_cluster_`, as.integer(family), as.integer(nb_cluster), theta, diffMax_NR, mu_in, lhs, sum_y, dum, obsCluster, table, cumtable, as.integer(nthreads))
}

get_n_cells <- function(index_i, index_j) {
  .Call(`_fixest2_get_n_cells_`, as.integer(index_i), as.integer(index_j))
}

cpp_derivconv_seq_gnl <- function(iterMax, diffMax, n_vars, nb_cluster_all, ll_d2, jacob_vector, deriv_init_vector, dum_vector) {
  .Call(`_fixest2_cpp_derivconv_seq_gnl_`, as.integer(iterMax), diffMax, as.integer(n_vars), nb_cluster_all, ll_d2, jacob_vector, deriv_init_vector, dum_vector)
}

cpp_derivconv_acc_gnl <- function(iterMax, diffMax, n_vars, nb_cluster_all, ll_d2, jacob_vector, deriv_init_vector, dum_vector) {
  .Call(`_fixest2_cpp_derivconv_acc_gnl_`, as.integer(iterMax), diffMax, as.integer(n_vars), nb_cluster_all, ll_d2, jacob_vector, deriv_init_vector, dum_vector)
}

cpp_derivconv_acc_2 <- function(iterMax, diffMax, n_vars, nb_cluster_all, n_cells, index_i, index_j, ll_d2, order, jacob_vector, deriv_init_vector, dum_vector) {
  .Call(`_fixest2_cpp_derivconv_acc_2_`, as.integer(iterMax), diffMax, as.integer(n_vars), nb_cluster_all, as.integer(n_cells), index_i, index_j, ll_d2, order, jacob_vector, deriv_init_vector, dum_vector)
}

cpp_derivconv_seq_2 <- function(iterMax, diffMax, n_vars, nb_cluster_all, n_cells, index_i, index_j, order, ll_d2, jacob_vector, deriv_init_vector, dum_vector) {
  .Call(`_fixest2_cpp_derivconv_seq_2_`, as.integer(iterMax), diffMax, as.integer(n_vars), nb_cluster_all, as.integer(n_cells), index_i, index_j, order, ll_d2, jacob_vector, deriv_init_vector, dum_vector)
}

update_deriv_single <- function(n_vars, nb_coef, r_ll_d2, r_jacob_vector, r_dum_vector) {
  .Call(`_fixest2_update_deriv_single_`, as.integer(n_vars), as.integer(nb_coef), r_ll_d2, r_jacob_vector, r_dum_vector)
}

cpp_conv_acc_gnl <- function(family, iterMax, diffMax, diffMax_NR, theta, nb_cluster_all, lhs, mu_init, dum_vector, tableCluster_vector, sum_y_vector, cumtable_vector, obsCluster_vector, nthreads = 1L) {
  .Call(`_fixest2_cpp_conv_acc_gnl_`, as.integer(family), as.integer(iterMax), diffMax, diffMax_NR, theta, nb_cluster_all, lhs, mu_init, dum_vector, tableCluster_vector, sum_y_vector, cumtable_vector, obsCluster_vector, as.integer(nthreads))
}

cpp_conv_seq_gnl <- function(family, iterMax, diffMax, diffMax_NR, theta, nb_cluster_all, lhs, mu_init, dum_vector, tableCluster_vector, sum_y_vector, cumtable_vector, obsCluster_vector, nthreads  = 1L) {
  .Call(`_fixest2_cpp_conv_seq_gnl_`, as.integer(family), as.integer(iterMax), diffMax, diffMax_NR, theta, nb_cluster_all, lhs, mu_init, dum_vector, tableCluster_vector, sum_y_vector, cumtable_vector, obsCluster_vector, as.integer(nthreads))
}

cpp_conv_acc_poi_2 <- function(n_i, n_j, n_cells, index_i, index_j, dum_vector, sum_y_vector, iterMax, diffMax, exp_mu_in, order) {
  .Call(`_fixest2_cpp_conv_acc_poi_2_`, as.integer(n_i), as.integer(n_j), as.integer(n_cells), index_i, index_j, dum_vector, sum_y_vector, as.integer(iterMax), diffMax, exp_mu_in, order)
}

cpp_conv_seq_poi_2 <- function(n_i, n_j, n_cells, index_i, index_j, dum_vector, sum_y_vector, iterMax, diffMax, exp_mu_in, order) {
  .Call(`_fixest2_cpp_conv_seq_poi_2_`, as.integer(n_i), as.integer(n_j), as.integer(n_cells), index_i, index_j, dum_vector, sum_y_vector, as.integer(iterMax), diffMax, exp_mu_in, order)
}

cpp_conv_acc_gau_2 <- function(n_i, n_j, n_cells, r_mat_row, r_mat_col, r_mat_value_Ab, r_mat_value_Ba, dum_vector, lhs, invTableCluster_vector, iterMax, diffMax, mu_in) {
  .Call(`_fixest2_cpp_conv_acc_gau_2_`, as.integer(n_i), as.integer(n_j), as.integer(n_cells), r_mat_row, r_mat_col, r_mat_value_Ab, r_mat_value_Ba, dum_vector, lhs, invTableCluster_vector, as.integer(iterMax), diffMax, mu_in)
}

cpp_conv_seq_gau_2 <- function(n_i, n_j, n_cells, r_mat_row, r_mat_col, r_mat_value_Ab, r_mat_value_Ba, dum_vector, lhs, invTableCluster_vector, iterMax, diffMax, mu_in) {
  .Call(`_fixest2_cpp_conv_seq_gau_2_`, as.integer(n_i), as.integer(n_j), as.integer(n_cells), r_mat_row, r_mat_col, r_mat_value_Ab, r_mat_value_Ba, dum_vector, lhs, invTableCluster_vector, as.integer(iterMax), diffMax, mu_in)
}

cpp_which_na_inf <- function(x, nthreads  = 1L) {
  .Call(`_fixest2_cpp_which_na_inf_`, x, as.integer(nthreads))
}

cpp_demean <- function(y, X_raw, r_weights, iterMax, diffMax, r_nb_id_Q, fe_id_list, table_id_I, slope_flag_Q, slope_vars_list, r_init, nthreads  = 1L, save_fixef = FALSE) {
  .Call(`_fixest2_cpp_demean_`, y, X_raw, r_weights, as.integer(iterMax), diffMax, r_nb_id_Q, fe_id_list, table_id_I, slope_flag_Q, slope_vars_list, r_init, as.integer(nthreads), save_fixef)
}

cpp_dsb <- function(Rstr) {
  .Call(`_fixest2_cpp_dsb_`, Rstr)
}

cpp_dsb_full_string <- function(Rstr) {
  .Call(`_fixest2_cpp_dsb_full_string_`, Rstr)
}

cpp_dsb_if_extract <- function(Rstr) {
  .Call(`_fixest2_cpp_dsb_if_extract_`, Rstr)
}

cpp_paste_conditional <- function(x, id, n) {
  .Call(`_fixest2_cpp_paste_conditional_`, x, as.integer(id), as.integer(n))
}

cpp_cholesky <- function(X, tol, nthreads  = 1L) {
  .Call(`_fixest2_cpp_cholesky_`, X, tol, as.integer(nthreads))
}

cpp_sparse_products <- function(X, w, y, correct_0w = FALSE, nthreads  = 1L) {
  .Call(`_fixest2_cpp_sparse_products_`, X, w, y, correct_0w, as.integer(nthreads))
}

cpppar_crossprod <- function(X, w, nthreads  = 1L) {
  .Call(`_fixest2_cpppar_crossprod_`, X, w, as.integer(nthreads))
}

cpp_mat_reconstruct <- function(X, id_excl) {
  .Call(`_fixest2_cpp_mat_reconstruct_`, X, id_excl)
}

cpp_iv_products <- function(X, y, Z, u, w, nthreads  = 1L) {
  .Call(`_fixest2_cpp_iv_products_`, X, y, Z, u, w, as.integer(nthreads))
}

cpp_iv_product_completion <- function(XtX, Xty, X, y, U, w, nthreads  = 1L) {
  .Call(`_fixest2_cpp_iv_product_completion_`, XtX, Xty, X, y, U, w, as.integer(nthreads))
}

cpp_iv_resid <- function(resid_2nd, coef, resid_1st, is_int, nthreads  = 1L) {
  .Call(`_fixest2_cpp_iv_resid_`, resid_2nd, coef, resid_1st, is_int, as.integer(nthreads))
}

cpp_add_commas <- function(x, r = 1L, whole = TRUE) {
  .Call(`_fixest2_cpp_add_commas_`, x, as.integer(r), whole)
}

cpp_find_never_always_treated <- function(cohort, period) {
  .Call(`_fixest2_cpp_find_never_always_treated_`, as.integer(cohort), period)
}

cpp_get_first_item <- function(x, n_items) {
  .Call(`_fixest2_cpp_get_first_item_`, as.integer(x), as.integer(n_items))
}

cpp_combine_clusters <- function(cluster_list, index) {
  .Call(`_fixest2_cpp_combine_clusters_`, cluster_list, as.integer(index))
}

cpp_cut <- function(x_sorted, cut_points, is_included) {
  .Call(`_fixest2_cpp_cut_`, x_sorted, cut_points, as.integer(is_included))
}

cpp_is_int <- function(x) {
  .Call(`_fixest2_cpp_is_int_`, x)
}

cpp_hash_string <- function(x) {
  .Call(`_fixest2_cpp_hash_string_`, x)
}

cpp_isConstant <- function(x) {
  .Call(`_fixest2_cpp_isConstant_`, x)
}

cpp_any_na_null <- function(x) {
  .Call(`_fixest2_cpp_any_na_null_`, x)
}

cpp_find_duplicates <- function(id, time) {
  .Call(`_fixest2_cpp_find_duplicates_`, as.integer(id), as.integer(time))
}

cpp_check_nested <- function(fe_list, cluster_list, fe_sizes, n) {
  .Call(`_fixest2_cpp_check_nested_`, fe_list, cluster_list, as.integer(fe_sizes), as.integer(n))
}

cpp_pgcd <- function(x) {
  .Call(`_fixest2_cpp_pgcd_`, as.integer(x))
}

cpp_get_fe_gnl <- function(Q, N, sumFE, dumMat, cluster_sizes, obsCluster) {
  storage.mode(dumMat) <- "integer"
  storage.mode(obsCluster) <- "integer"
  .Call(`_fixest2_cpp_get_fe_gnl_`, as.integer(Q), as.integer(N), sumFE, dumMat, as.integer(cluster_sizes), obsCluster)
}

cpp_escape_markup <- function(Rstr) {
  .Call(`_fixest2_cpp_escape_markup_`, Rstr)
}

cpp_factor_matrix <- function(fact, is_na_all, who_is_dropped, var, col_names) {
  .Call(`_fixest2_cpp_factor_matrix_`, as.integer(fact), is_na_all, as.integer(who_is_dropped), var, col_names)
}

cpp_diag_XUtX <- function(X, U) {
  .Call(`_fixest2_cpp_diag_XUtX_`, X, U)
}

cpp_lgamma <- function(x) {
  .Call(`_fixest2_cpp_lgamma_`, x)
}

cpp_log_a_exp <- function(a, mu, exp_mu) {
  .Call(`_fixest2_cpp_log_a_exp_`, a, mu, exp_mu)
}

cpp_partialDerivative_other <- function(iterMax, Q, N, epsDeriv, ll_d2, dx_dother, init, dumMat, nbCluster) {
  .Call(`_fixest2_cpp_partialDerivative_other_`, as.integer(iterMax), as.integer(Q), as.integer(N), epsDeriv, ll_d2, dx_dother, init, dumMat, as.integer(nbCluster))
}

cpp_table <- function(Q, dum) {
  .Call(`_fixest2_cpp_table_`, as.integer(Q), as.integer(dum))
}

cpp_ssr_null <- function(y, w = 0) {
  .Call(`_fixest2_cpp_ssr_null_`, y, w)
}

cpp_ssq <- function(x, w = 0) {
  .Call(`_fixest2_cpp_ssq_`, x, w)
}

cpp_constant_dum <- function(k, x, dum, only_0) {
  .Call(`_fixest2_cpp_constant_dum_`, as.integer(k), x, as.integer(dum), only_0 = FALSE)
}

cpp_tapply_sum <- function(Q, x, dum) {
  .Call(`_fixest2_cpp_tapply_sum_`, as.integer(Q), x, as.integer(dum))
}

cpp_tapply_vsum <- function(Q, x, dum) {
  .Call(`_fixest2_cpp_tapply_vsum_`, as.integer(Q), x, as.integer(dum))
}

cpp_lag_obs <- function(id, time, nlag) {
  .Call(`_fixest2_cpp_lag_obs_`, as.integer(id), as.integer(time), as.integer(nlag))
}

cpp_quf_gnl <- function(x) {
  .Call(`_fixest2_cpp_quf_gnl_`, x)
}

cpppar_quf_table_sum <- function(x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, nthreads = 1L, do_refactor, r_x_sizes, obs2keep) {
  .Call(`_fixest2_cpppar_quf_table_sum_`, x, y, do_sum_y, rm_0, rm_1, rm_single, only_slope, as.integer(nthreads), do_refactor, r_x_sizes, as.integer(obs2keep))
}

cpp_get_nb_threads <- function() {
  .Call(`_fixest2_cpp_get_nb_threads`)
}

cpppar_which_na_inf_vec <- function(x, nthreads = 1L) {
  .Call(`_fixest2_cpppar_which_na_inf_vec_`, x, as.integer(nthreads))
}

cpppar_which_na_inf_mat <- function(mat, nthreads = 1L) {
  .Call(`_fixest2_cpppar_which_na_inf_mat_`, mat, as.integer(nthreads))
}

cpppar_which_na_inf_df <- function(df, nthreads = 1L) {
  .Call(`_fixest2_cpppar_which_na_inf_df_`, df, as.integer(nthreads))
}

cpppar_check_only_0 <- function(x_mat, nthreads = 1L) {
  .Call(`_fixest2_cpppar_check_only_0_`, x_mat, as.integer(nthreads))
}

cpppar_exp <- function(x, nthreads = 1L) {
  .Call(`_fixest2_cpppar_exp_`, x, as.integer(nthreads))
}

cpppar_log <- function(x, nthreads = 1L) {
  .Call(`_fixest2_cpppar_log_`, x, as.integer(nthreads))
}

cpppar_log_a_exp <- function(nthreads = 1L, a, mu, exp_mu) {
  .Call(`_fixest2_cpppar_log_a_exp_`, as.integer(as.integer(nthreads)), a, mu, exp_mu)
}

cpppar_cond_means <- function(mat_vars, treat, nthreads = 1L) {
  .Call(`_fixest2_cpppar_cond_means_`, mat_vars, as.integer(treat), as.integer(nthreads))
}

cpppar_lgamma <- function(x, nthreads = 1L) {
  .Call(`_fixest2_cpppar_lgamma_`, x, as.integer(nthreads))
}

cpppar_digamma <- function(x, nthreads = 1L) {
  .Call(`_fixest2_cpppar_digamma_`, x, as.integer(nthreads))
}

cpppar_trigamma <- function(x, nthreads = 1L) {
  .Call(`_fixest2_cpppar_trigamma_`, x, as.integer(nthreads))
}

cpppar_poisson_linkinv <- function(x, nthreads = 1L) {
  .Call(`_fixest2_cpppar_poisson_linkinv_`, x, as.integer(nthreads))
}

cpppar_poisson_validmu <- function(x, nthreads = 1L) {
  .Call(`_fixest2_cpppar_poisson_validmu_`, x, as.integer(nthreads))
}

cpppar_logit_linkfun <- function(x, nthreads = 1L) {
  .Call(`_fixest2_cpppar_logit_linkfun_`, x, as.integer(nthreads))
}

cpppar_logit_linkinv <- function(x, nthreads = 1L) {
  .Call(`_fixest2_cpppar_logit_linkinv_`, x, as.integer(nthreads))
}

cpppar_logit_mueta <- function(x, nthreads = 1L) {
  .Call(`_fixest2_cpppar_logit_mueta_`, x, as.integer(nthreads))
}

cpppar_logit_devresids <- function(y, mu, wt, nthreads = 1L) {
  .Call(`_fixest2_cpppar_logit_devresids_`, y, mu, wt, as.integer(nthreads))
}

cpppar_xwy <- function(X, y, w, nthreads = 1L) {
  .Call(`_fixest2_cpppar_xwy_`, X, y, w, as.integer(nthreads))
}

cpppar_xbeta <- function(X, beta, nthreads = 1L) {
  .Call(`_fixest2_cpppar_xbeta_`, X, beta, as.integer(nthreads))
}

cpppar_matprod <- function(x, y, nthreads = 1L) {
  .Call(`_fixest2_cpppar_matprod_`, x, y, as.integer(nthreads))
}

cpp_colon_to_star <- function(Rstr) {
  .Call(`_fixest2_cpp_colon_to_star_`, Rstr)
}

cpp_newey_west <- function(S, w, nthreads = 1L) {
  .Call(`_fixest2_cpp_newey_west_`, S, w, as.integer(nthreads))
}

cpp_newey_west_panel <- function(S, w, unit, G, time, T, nthreads = 1L) {
  .Call(`_fixest2_cpp_newey_west_panel_`, S, w, as.integer(unit), as.integer(G), as.integer(time), as.integer(T), as.integer(nthreads))
}

cpp_driscoll_kraay <- function(S, w, time, T, nthreads = 1L) {
  .Call(`_fixest2_cpp_driscoll_kraay_`, S, w, as.integer(time), as.integer(T), as.integer(nthreads))
}

cpp_vcov_conley <- function(S, lon_rad, lat_rad, distance, cutoff, nthreads = 1L) {
  .Call(`_fixest2_cpp_vcov_conley_`, S, lon_rad, lat_rad, as.integer(distance), cutoff, as.integer(nthreads))
}
