// Demeans each variable in input
// The method is based on obtaining the optimal cluster coefficients
// Works but only the master thread can call that function!
// What happens if the master thread has finished its job but the lower thread
// is in an "infinite" loop? this is tricky, as such we cannot stop it.
// Solution: create a function keeping the threads idle waiting for the complete
// job to be done BUT I need to add static allocation of threads => performance
// cost

// List of objects, used to
// lighten the writting of the functions

#include "02_0_demeaning.hpp"

void FEClass::add_2_fe_coef_to_mu(double *fe_coef_a, double *fe_coef_b,
                                  double *in_out_C, double *out_N,
                                  bool update_beta = true) {
  // We add the value of the FE coefficients to each observation

  // Step 1: we update the coefficients of b

  if (update_beta) {
    compute_fe_coef_2_internal(fe_coef_a, fe_coef_b, in_out_C, out_N);
  }

  // Step 2: we add the value of each coef

  for (int q = 0; q < 2; ++q) {
    double *my_fe_coef = q == 0 ? fe_coef_a : fe_coef_b;
    int *my_fe = p_fe_id[q];
    bool is_slope = is_slope_Q[q];
    int V = nb_vs_Q[q];

    simple_mat_of_vs_vars VS_mat(this, q);
    simple_mat_with_id my_vs_coef(my_fe_coef, nb_vs_Q[q], 1);

    for (int i = 0; i < n_obs; ++i) {
      if (is_slope) {
        for (int v = 0; v < V; ++v) {
          out_N[i] += my_vs_coef(my_fe[i] - 1, v) * VS_mat(i, v);
        }
      } else {
        out_N[i] += my_fe_coef[my_fe[i] - 1];
      }
    }
  }
}

void demean_single_1(int v, PARAM_DEMEAN *args) {
  // v: variable identifier to demean

  // Q == 1: nothing to say, just compute the closed form

  // loading the data
  int nb_coef_T = args->nb_coef_T;

  vector<sVec> &p_input = args->p_input;
  vector<double *> &p_output = args->p_output;

  // fe_info
  FEClass &FE_info = *(args->p_FE_info);

  // vector of fixed-effects coefficients initialized at 0
  vector<double> fe_coef(nb_coef_T, 0);
  double *p_fe_coef = fe_coef.data();

  // interruption handling
  // bool isMaster = omp_get_thread_num() == 0;
  // bool *pStopNow = args->stopnow;
  // if (isMaster) {
  //   if (pending_interrupt()) {
  //     *pStopNow = true;
  //   }
  // }

  // the input & output
  sVec &input = p_input[v];
  double *output = p_output[v];

  // We compute the FEs
  FE_info.compute_fe_coef(p_fe_coef, input);

  // Output:
  FE_info.add_fe_coef_to_mu(0, p_fe_coef, output);

  // saving the fixef coefs
  double *fixef_values = args->fixef_values;
  if (args->save_fixef) {
    for (int m = 0; m < nb_coef_T; ++m) {
      fixef_values[m] = fe_coef[m];
    }
  }
}

bool demean_acc_gnl(int v, int iterMax, PARAM_DEMEAN *args,
                    bool two_fe = false) {
  //
  // data
  //

  // fe_info
  FEClass &FE_info = *(args->p_FE_info);

  // algo info
  const int n_extraProj = args->algo_extraProj;
  const int iter_projAfterAcc = args->algo_iter_projAfterAcc;
  const int iter_grandAcc = args->algo_iter_grandAcc;

  int n_obs = args->n_obs;
  int nb_coef_T = args->nb_coef_T;
  int Q = args->Q;
  double diffMax = args->diffMax;

  int nb_coef_all = nb_coef_T;
  const bool two_fe_algo = two_fe || Q == 2;
  if (two_fe_algo) {
    Q = 2;
    // special case: 2 FEs, below will contain only the size of the first FE
    //   this is because 2-FE is a special case which allows to avoid creating
    //   a N-size vector
    nb_coef_T = FE_info.nb_coef_Q[0];
    // nb_coef_all: only used in in-out
    nb_coef_all = FE_info.nb_coef_Q[0] + FE_info.nb_coef_Q[1];
  }

  // input output
  vector<sVec> &p_input = args->p_input;
  vector<double *> &p_output = args->p_output;
  sVec &input = p_input[v];
  double *output = p_output[v];

  // temp var:
  int size_other_means = two_fe_algo ? FE_info.nb_coef_Q[1] : n_obs;
  vector<double> sum_other_means_or_second_coef(size_other_means);
  double *p_sum_other_means = sum_other_means_or_second_coef.data();

  // conditional sum of input minus output
  vector<double> sum_input_output(nb_coef_all, 0);
  double *p_sum_in_out = sum_input_output.data();

  for (int q = 0; q < Q; ++q) {
    FE_info.compute_in_out(q, p_sum_in_out, input, output);
  }

  // interruption handling
  // bool isMaster = omp_get_thread_num() == 0;
  // bool *pStopNow = args->stopnow;
  // I overcast to remember the lesson
  // rough estimate nber operation per iter
  // double flop = 4.0 * (5 + 12 * (Q - 1) + 4 * (Q - 1) * (Q - 1)) *
  //               static_cast<double>(n_obs);
  // if (two_fe_algo) {
  //   flop = 20.0 * static_cast<double>(n_obs);
  // }
  // int iterSecond =
  //     ceil(2000000000 / flop / 5);  // nber iter per 1/5 second at 2GHz

  //
  // IT iteration (preparation)
  //

  // variables on 1:nb_coef
  // note that in the case of two_fe_algo, the length is equal to the 1st FE
  // coefs
  vector<double> X(nb_coef_T, 0);
  vector<double> GX(nb_coef_T);
  vector<double> GGX(nb_coef_T);
  // pointers:
  double *p_X = X.data();
  double *p_GX = GX.data();
  double *p_GGX = GGX.data();

  // variables on 1:(Q-1)
  int nb_coef_no_Q = 0;
  for (int q = 0; q < (Q - 1); ++q) {
    nb_coef_no_Q += FE_info.nb_coef_Q[q];
  }
  vector<double> delta_GX(nb_coef_no_Q);
  vector<double> delta2_X(nb_coef_no_Q);

  // additional vectors used to save the coefs
  vector<double> Y(nb_coef_T);
  vector<double> GY(nb_coef_T);
  vector<double> GGY(nb_coef_T);
  // pointers:
  double *p_Y = Y.data();
  double *p_GY = GY.data();
  double *p_GGY = GGY.data();

  int grand_acc = 0;

  //
  // the main loop
  //

  // first iteration
  compute_fe(Q, p_X, p_GX, p_sum_other_means, p_sum_in_out, args);

  // check whether we should go into the loop
  bool keepGoing = false;
  for (int i = 0; i < nb_coef_T; ++i) {
    if (continue_crit(X[i], GX[i], diffMax)) {
      keepGoing = true;
      break;
    }
  }

  // For the stopping criterion on total addition
  double ssr = 0;

  int iter = 0;
  bool numconv = false;
  // while (!*pStopNow && keepGoing && iter < iterMax) {
  while (keepGoing && iter < iterMax) {
    // if (isMaster && iter % iterSecond == 0) {
    //   if (pending_interrupt()) {
    //     *pStopNow = true;
    //     break;
    //   }
    // }

    iter++;

    for (int rep = 0; rep < n_extraProj; ++rep) {
      // simple projections, at the request of the user
      // may be useful to hasten convergence on special cases
      // default is 0
      compute_fe(Q, p_GX, p_GGX, p_sum_other_means, p_sum_in_out, args);
      compute_fe(Q, p_GGX, p_X, p_sum_other_means, p_sum_in_out, args);
      compute_fe(Q, p_X, p_GX, p_sum_other_means, p_sum_in_out, args);
    }

    // origin: GX, destination: GGX
    compute_fe(Q, p_GX, p_GGX, p_sum_other_means, p_sum_in_out, args);

    // X: outcome of the acceleration
    numconv =
        dm_update_X_IronsTuck(nb_coef_no_Q, X, GX, GGX, delta_GX, delta2_X);
    if (numconv) break;

    if (iter >= iter_projAfterAcc) {
      memcpy(p_Y, p_X, nb_coef_T * sizeof(double));
      compute_fe(Q, p_Y, p_X, p_sum_other_means, p_sum_in_out, args);
    }

    // origin: X, destination: GX
    compute_fe(Q, p_X, p_GX, p_sum_other_means, p_sum_in_out, args);

    keepGoing = false;
    for (int i = 0; i < nb_coef_no_Q; ++i) {
      if (continue_crit(X[i], GX[i], diffMax)) {
        keepGoing = true;
        break;
      }
    }

    if (iter % iter_grandAcc == 0) {
      ++grand_acc;
      if (grand_acc == 1) {
        memcpy(p_Y, p_GX, nb_coef_T * sizeof(double));
      } else if (grand_acc == 2) {
        memcpy(p_GY, p_GX, nb_coef_T * sizeof(double));
      } else {
        memcpy(p_GGY, p_GX, nb_coef_T * sizeof(double));
        numconv =
            dm_update_X_IronsTuck(nb_coef_no_Q, Y, GY, GGY, delta_GX, delta2_X);
        if (numconv) break;
        compute_fe(Q, p_Y, p_GX, p_sum_other_means, p_sum_in_out, args);
        grand_acc = 0;
      }
    }

    // Other stopping criterion: change to SSR very small
    if (iter % 40 == 0) {
      // mu_current is the vector of means
      vector<double> mu_current(n_obs, 0);
      double *p_mu = mu_current.data();

      if (two_fe_algo) {
        FE_info.add_2_fe_coef_to_mu(p_GX, p_sum_other_means, p_sum_in_out,
                                    p_mu);
      } else {
        for (int q = 0; q < Q; ++q) {
          FE_info.add_fe_coef_to_mu(q, p_GX, p_mu);
        }
      }

      double ssr_old = ssr;

      // we compute the new SSR
      ssr = 0;
      double resid;
      for (int i = 0; i < n_obs; ++i) {
        resid = input[i] - mu_current[i];
        ssr += resid * resid;
      }

      if (stopping_crit(ssr_old, ssr, diffMax)) {
        break;
      }
    }
  }

  //
  // Updating the output
  //

  double *p_beta_final = nullptr;
  if (two_fe_algo) {
    // we end with a last iteration
    double *p_alpha_final = p_GX;
    p_beta_final = p_sum_other_means;

    FE_info.compute_fe_coef_2(p_alpha_final, p_alpha_final, p_beta_final,
                              p_sum_in_out);
    FE_info.add_2_fe_coef_to_mu(p_alpha_final, p_beta_final, p_sum_in_out,
                                output, false);
  } else {
    for (int q = 0; q < Q; ++q) {
      FE_info.add_fe_coef_to_mu(q, p_GX, output);
    }
  }

  // keeping track of iterations
  int *iterations_all = args->p_iterations_all;
  iterations_all[v] += iter;

  // saving the fixef coefs
  double *fixef_values = args->fixef_values;
  if (args->save_fixef) {
    for (int m = 0; m < nb_coef_T; ++m) {
      fixef_values[m] += GX[m];
    }

    if (two_fe_algo) {
      // we add the other coefficients
      int n_coefs_FE1 = nb_coef_T;
      int n_coefs_FE2 = size_other_means;
      for (int m = 0; m < n_coefs_FE2; ++m) {
        fixef_values[n_coefs_FE1 + m] += p_beta_final[m];
      }
    }
  }

  bool conv = iter == iterMax ? false : true;

  return (conv);
}

void demean_single_gnl(int v, PARAM_DEMEAN *args) {
  // v: integer identifying the variable to demean

  // Algorithm to quickly get the means
  // Q >= 3 => acceleration for 15 iter
  // if no convergence: conv 2 FEs
  // then acceleration again

  // data
  int iterMax = args->iterMax;

  // Note that the historical default was 15 iterations (still the case on
  // 2024-02-14)
  int iter_warmup = args->algo_iter_warmup;
  int Q = args->Q;

  if (Q == 2) {
    demean_acc_gnl(v, iterMax, args);
  } else {
    bool conv = false;

    // iter_warmup <= 0 means no warmup and we start directly with 2FE algo
    if (iter_warmup > 0) {
      conv = demean_acc_gnl(v, iter_warmup, args);
      iter_warmup = 0;
    }

    if (conv == false && iterMax > iter_warmup) {
      // convergence for the first 2 FEs
      int iter_max_2FE = iterMax / 2 - iter_warmup;
      if (iter_max_2FE > 0) {
        demean_acc_gnl(v, iter_max_2FE, args, true);
      }

      // re-acceleration
      int iter_previous = args->p_iterations_all[v];
      demean_acc_gnl(v, iterMax - iter_previous, args);
    }
  }

  int *jobdone = args->jobdone;
  jobdone[v] = 1;
}

// Loop over demean_single
[[cpp11::register]] list cpp_demean_(
    SEXP y, SEXP X_raw, SEXP r_weights, int iterMax, double diffMax,
    SEXP r_nb_id_Q, SEXP fe_id_list, SEXP table_id_I, SEXP slope_flag_Q,
    SEXP slope_vars_list, SEXP r_init, int nthreads, int algo_extraProj,
    int algo_iter_warmup, int algo_iter_projAfterAcc, int algo_iter_grandAcc,
    bool save_fixef) {
  // main fun that calls demean_single
  // preformats all the information needed on the fixed-effects
  // y: the dependent variable
  // X_raw: the matrix of the explanatory variables -- can be "empty"

  // when including weights: recreate table values
  // export weights and is_weight bool

  // slope_flag: whether a FE is a varying slope
  // slope_var: the associated variables with varying slopes

  // initial variables
  int Q = Rf_length(r_nb_id_Q);

  // info on y
  sMat m_y(y);
  int n_vars_y = m_y.ncol();
  bool useY = n_vars_y > 0;
  bool is_y_list = n_vars_y > 1 || TYPEOF(y) == VECSXP;
  int n_obs = m_y.nrow();

  // info on X
  sMat m_X(X_raw);
  int n_vars_X = m_X.ncol();
  if (useY == false) {
    n_obs = m_X.nrow();
  }
  bool useX = n_vars_X > 0;

  if (n_obs == 0 || n_obs == 1) {
    // The data set is of length 1!!!!
    n_obs = 1;

    m_y = sMat(y, true);
    n_vars_y = m_y.ncol();
    useY = true;

    m_X = sMat(X_raw, true);
    n_vars_X = m_X.ncol();
    useX = n_vars_X > 0;
  }

  int n_vars = n_vars_y + n_vars_X;

  // initialisation if needed (we never initialize when only one FE, except if
  // asked explicitly)
  bool isInit = Rf_xlength(r_init) != 1 && Q > 1;
  double *init = REAL(r_init);
  bool saveInit = ((isInit || init[0] != 0) && Q > 1) || init[0] == 666;

  // Creating the object containing all information on the FEs
  FEClass FE_info(n_obs, Q, r_weights, fe_id_list, r_nb_id_Q, table_id_I,
                  slope_flag_Q, slope_vars_list);
  int nb_coef_T = FE_info.nb_coef_T;

  // output vector: (Note that if the means are provided, we use that vector and
  // will modify it in place)
  int64_t n_total = static_cast<int64_t>(n_obs) * n_vars;
  SEXP output_values = PROTECT(Rf_allocVector(REALSXP, isInit ? 1 : n_total));

  double *p_output_origin = isInit ? init : REAL(output_values);
  if (isInit == false) {
    std::fill_n(p_output_origin, n_total, 0);
  }

  // vector of pointers: input/output

  vector<double *> p_output(n_vars);
  p_output[0] = p_output_origin;
  for (int v = 1; v < n_vars; v++) {
    p_output[v] = p_output[v - 1] + n_obs;
  }

  vector<sVec> p_input(n_vars);

  for (int k = 0; k < n_vars_X; ++k) {
    p_input[k] = m_X[k];
  }

  for (int k = 0; k < n_vars_y; ++k) {
    p_input[n_vars_X + k] = m_y[k];
  }

  // keeping track of iterations
  vector<int> iterations_all(n_vars, 0);
  int *p_iterations_all = iterations_all.data();

  // save fixef option
  if (useX && save_fixef) {
    stop("save_fixef can only be used when there is no Xs.");
  }

  vector<double> fixef_values(save_fixef ? nb_coef_T : 1, 0);
  double *p_fixef_values = fixef_values.data();

  // Sending variables to envir

  PARAM_DEMEAN args;

  args.n_obs = n_obs;
  args.iterMax = iterMax;
  args.diffMax = diffMax;
  args.Q = Q;
  args.nb_coef_T = nb_coef_T;
  args.p_input = p_input;
  args.p_output = p_output;
  args.p_iterations_all = p_iterations_all;
  args.algo_extraProj = algo_extraProj;
  args.algo_iter_warmup = algo_iter_warmup;
  args.algo_iter_projAfterAcc = algo_iter_projAfterAcc;
  // negative number or 0 means never
  args.algo_iter_grandAcc =
      algo_iter_grandAcc <= 0 ? 1000000 : algo_iter_grandAcc;

  // save FE_info
  args.p_FE_info = &FE_info;

  // save fixef:
  args.save_fixef = save_fixef;
  args.fixef_values = p_fixef_values;

  // stopping flag + indicator that job is finished
  bool stopnow = false;
  args.stopnow = &stopnow;
  vector<int> jobdone(n_vars, 0);
  int *pjobdone = jobdone.data();
  args.jobdone = pjobdone;

  int counter = 0;
  int *pcounter = &counter;

  // main loop

  int nthreads_current = nthreads > n_vars ? n_vars : nthreads;

#pragma omp parallel for num_threads(nthreads_current) schedule(static, 1)
  for (int v = 0; v < (n_vars + nthreads_current); ++v) {
    // demean_single is the workhorse
    // you get the "mean"

    if (!*(args.stopnow)) {
      if (v < n_vars) {
        if (Q == 1) {
          demean_single_1(v, &args);
        } else {
          demean_single_gnl(v, &args);
        }
      } else if (true && Q != 1) {
        stayIdleCheckingInterrupt(&stopnow, jobdone, n_vars, pcounter);
      }
    }
  }

  if (*(args.stopnow)) {
    UNPROTECT(1);
    stop("cpp_demean: User interrupt.");
  }

  // save

  writable::list res;  // a vector and a matrix

  int nrow = useX ? n_obs : 1;
  int ncol = useX ? n_vars_X : 1;
  writable::doubles_matrix<> X_demean(nrow, ncol);
  // TODO: check if I really need to fill with 0s
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      X_demean(i, j) = 0;
    }
  }

  sVec p_input_tmp;
  double *p_output_tmp;
  for (int k = 0; k < n_vars_X; ++k) {
    p_input_tmp = p_input[k];
    p_output_tmp = p_output[k];

    for (int i = 0; i < n_obs; ++i) {
      X_demean(i, k) = p_input_tmp[i] - p_output_tmp[i];
    }
  }

  res.push_back({"X_demean"_nm = X_demean});

  if (is_y_list && useY) {
    writable::list y_demean(n_vars_y);

    for (int v = 0; v < n_vars_y; ++v) {
      p_input_tmp = p_input[n_vars_X + v];
      p_output_tmp = p_output[n_vars_X + v];

      writable::doubles y_demean_tmp(n_obs);
      for (int i = 0; i < n_obs; ++i) {
        y_demean_tmp[i] = p_input_tmp[i] - p_output_tmp[i];
      }

      y_demean[v] = y_demean_tmp;
    }

    res.push_back({"y_demean"_nm = y_demean});
  } else {
    writable::doubles y_demean(useY ? n_obs : 1);
    if (useY) {
      // y is always the last variable
      p_input_tmp = p_input[n_vars - 1];
      p_output_tmp = p_output[n_vars - 1];
      for (int i = 0; i < n_obs; ++i) {
        y_demean[i] = p_input_tmp[i] - p_output_tmp[i];
      }
    }
    res.push_back({"y_demean"_nm = y_demean});
  }

  // iterations
  writable::integers iter_final(n_vars);
  for (int v = 0; v < n_vars; ++v) {
    iter_final[v] = p_iterations_all[v];
  }

  res.push_back({"iterations"_nm = iter_final});

  // if save is requested
  if (saveInit) {
    if (isInit) {
      res.push_back({"means"_nm = r_init});
    } else {
      res.push_back({"means"_nm = output_values});
    }
  } else {
    res.push_back({"means"_nm = 0.0});
  }

  // save fixef coef
  if (save_fixef) {
    res.push_back({"fixef_coef"_nm = fixef_values});
  }

  UNPROTECT(1);

  return (res);
}
