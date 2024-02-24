#include "06_0_quf.hpp"

void quf_single(void *px_in, std::string &x_type, int n, int *x_uf,
                std::vector<double> &x_unik) {
  // preparation for strings
  bool IS_STR = x_type == "string";

  bool IS_INT = false;
  bool is_int_in_double = false;
  if (x_type == "double") {
    // we check if underlying structure is int
    IS_INT = true;
    double *px = (double *)px_in;
    for (int i = 0; i < n; ++i) {
      if (!(px[i] == (int)px[i])) {
        IS_INT = false;
        break;
      }
    }

    is_int_in_double = IS_INT;  // true: only if x is REAL + INT test OK

  } else if (x_type == "int") {
    IS_INT = true;
  }

  if (IS_INT) {
    // integer

    double max_value;
    void *px_generic;
    int X_MIN;
    if (is_int_in_double) {
      double *px = (double *)px_in;
      px_generic = (double *)px_in;
      double x_min = px[0], x_max = px[0], x_tmp;
      for (int i = 1; i < n; ++i) {
        x_tmp = px[i];
        if (x_tmp > x_max) x_max = x_tmp;
        if (x_tmp < x_min) x_min = x_tmp;
      }
      X_MIN = static_cast<int>(x_min);
      max_value = x_max - x_min;
    } else {
      int *px = (int *)px_in;
      px_generic = (int *)px_in;
      int x_min = px[0], x_max = px[0], x_tmp;
      for (int i = 1; i < n; ++i) {
        x_tmp = px[i];
        if (x_tmp > x_max) x_max = x_tmp;
        if (x_tmp < x_min) x_min = x_tmp;
      }
      X_MIN = x_min;
      max_value = x_max - x_min;
    }

    // creating + copying a 10**5 vector takes about 0.5ms which is OK even if n
    // << 10**5 quf_int is really extremely efficient, so we use it whenever
    // appropriate In quf_int_gnl we create two n-size vectors + the method is
    // less efficient by construction so we go into quf_int whenever max_value
    // <= 2.5*n

    if (max_value < 100000 || max_value <= 2.5 * n) {
      quf_int(n, x_uf, px_generic, x_unik, X_MIN, static_cast<int>(max_value),
              is_int_in_double);
    } else if (max_value < 0x10000000) {
      // we don't cover ranges > 2**31 (uints are pain in the neck)
      quf_int_gnl(n, x_uf, px_generic, x_unik, X_MIN, is_int_in_double);
    } else {
      // ranges > 2**31 => as double

      if (is_int_in_double) {
        quf_double(n, x_uf, (double *)px_generic, x_unik, false);
      } else {
        // we need to create a vector of double, otherwise: pointer issue
        std::vector<double> x_dble(n);
        int *px = (int *)px_in;
        for (int i = 0; i < n; ++i) x_dble[i] = static_cast<double>(px[i]);
        quf_double(n, x_uf, x_dble.data(), x_unik, false);
      }
    }

  } else if (IS_STR) {
    // string -- beforehand transformed as ULL
    quf_double(n, x_uf, px_in, x_unik, true);
  } else {
    // double
    quf_double(n, x_uf, px_in, x_unik, false);
  }
}

void quf_refactor(int *px_in, int x_size, integers &obs2keep, int n, int *x_uf,
                  std::vector<double> &x_unik, std::vector<int> &x_table) {
  // pxin: data that has been qufed already => so integer ranging from 1 to nber
  // of unique elements obs2keep => optional, observations to keep

  int n_keep = obs2keep.size();
  bool keep_obs = obs2keep[0] != 0;

  if (keep_obs) {
    std::vector<int> id_new(x_size, 0);
    int val = 0;
    int val_new = 1;
    for (int i = 0; i < n_keep; ++i) {
      val = px_in[obs2keep[i] - 1] - 1;
      if (id_new[val] == 0) {
        x_table.push_back(1);
        x_unik.push_back(val + 1);
        id_new[val] = val_new++;
      } else {
        ++x_table[id_new[val] - 1];
      }
      x_uf[i] = id_new[val];
    }

  } else {
    // We just create the table and the unique
    x_table.resize(x_size);
    std::fill(x_table.begin(), x_table.end(), 0);

    for (int i = 0; i < n; ++i) {
      ++x_table[px_in[i] - 1];
    }

    x_unik.resize(x_size);
    for (int i = 0; i < x_size; ++i) {
      x_unik[i] = i + 1;
    }
  }
}

void quf_refactor_table_sum_single(
    int n, int *quf_old, int *quf_new, std::vector<bool> &obs_removed,
    std::vector<double> &x_unik, std::vector<double> &x_unik_new,
    std::vector<double> &x_removed, std::vector<int> &x_table, double *py,
    std::vector<double> &sum_y, bool do_sum_y, bool rm_1, bool rm_single,
    std::vector<bool> &id_pblm, bool check_pblm, bool *pstop_now, int n_FE) {
  // takes in the old quf, the observations removed,
  // then recreates the vectors of:
  // quf, table, sum_y, x_unik_new

  // IDs range from 1 to n_values

  int D = x_unik.size();

  // Only if rm_1 (0s and 1s removed) that we recreate the pblmatic IDs
  // because if only 0s, only the ones in the original id_pblm are removed
  // we find out which ones are only 0s or only 1s
  // NOTA: I don't reiterate the observation removal => just one round
  //  (otherwise it could go several rounds => too time consuming, no big value
  //  added)

  // if n_FE = 1, no need to recheck the removal since by definition it can't be
  // messed up by the other dimensions. Everything is done properly directly.

  if ((rm_single || rm_1 || !check_pblm) && !*pstop_now && n_FE > 1) {
    // If we haven't excluded observations because of 0 only values (ie
    // check_pblm == false) then we here need to update id_pblm, because some
    // IDs could have been taken out due to other FEs being removed we also need
    // to recheck when singletons are removed

    // The fact that we don't need to recheck when removing only 0s is an oddity
    // Indeed, when removing only 0s, the IDs of one dimension aren't messed up
    // by the removal of observations from other dimensions because if they
    // were to be removed because of this, this would mean that all values
    // of this ID are 0, so they'd be removed anyway from their own dimension
    // in the first place.

    // the new ID problem => what are the IDs that still exists?
    std::vector<bool> id_still_exists(D, false);

    // We also check that some observations are still left!

    // we recompute the table and sum_y
    std::fill(x_table.begin(), x_table.end(), 0);
    // std::fill(sum_y.begin(), sum_y.end(), 0);

    if (do_sum_y) {
      std::fill(sum_y.begin(), sum_y.end(), 0);
    }

    int id = 0;
    for (int i = 0; i < n; ++i) {
      if (!obs_removed[i]) {
        id = quf_old[i] - 1;
        if (!id_still_exists[id]) id_still_exists[id] = true;
        ++x_table[id];
        if (do_sum_y) sum_y[id] += py[i];
      }
    }

    // recreating id_pblm
    id_pblm.resize(D);  // we resize because we could have had length = 0
    for (int d = 0; d < D; ++d) id_pblm[d] = !id_still_exists[d];

    // Loop finding out the problems
    std::vector<bool> id_pblm_check(D, false);
    if (rm_1) {
      for (int d = 0; d < D; ++d) {
        if (sum_y[d] == 0 || sum_y[d] == x_table[d]) {
          id_pblm_check[d] = true;
        }
      }
    }

    int sum_problems =
        std::accumulate(id_pblm_check.begin(), id_pblm_check.end(), 0);
    // Rcout << "sum_problems = " << sum_problems << "\n";
    *pstop_now = sum_problems == D;
  }

  if (*pstop_now == false) {
    // is there a problem?
    bool any_pblm = false;
    // we check only if needed (if !rm_single && !rm_1 && length(id_pblm) == 0
    // => means no problem)
    if (!id_pblm.empty()) {
      for (int d = 0; d < D; ++d) {
        if (id_pblm[d]) {
          any_pblm = true;
          break;
        }
      }
    }

    // Creating the new IDs
    std::vector<int> id_new;
    int nb_pblm = 0;
    if (any_pblm) {
      // we create these new IDs only if there is a pblm
      // if no problem then the elements of quf_new are the same as in quf_old

      id_new.resize(D);
      std::fill(id_new.begin(), id_new.end(), 0);

      for (int d = 0; d < D; ++d) {
        if (id_pblm[d]) {
          ++nb_pblm;
        } else {
          id_new[d] = 1 + d - nb_pblm;
        }
      }
    }

    // New quf + table + sum_y

    int D_new = D - nb_pblm;

    // Rcout << "D new: " << D_new << "\n";
    // Rcout << "any_pblm: " << any_pblm << "\n";
    // Rcout << "id_pblm size: " << id_pblm.size() << "\n";

    x_table.resize(D_new);
    std::fill(x_table.begin(), x_table.end(), 0);
    if (do_sum_y) {
      sum_y.resize(D_new);
      std::fill(sum_y.begin(), sum_y.end(), 0);
    }

    int id;
    int i_new = 0;
    for (int i = 0; i < n; ++i) {
      if (!obs_removed[i]) {
        id = any_pblm ? id_new[quf_old[i] - 1] : quf_old[i];
        ++x_table[id - 1];
        if (do_sum_y) sum_y[id - 1] += py[i];
        quf_new[i_new++] = id;
      }
    }

    // x_unik + x_removed

    if (any_pblm) {
      // we create x_unik_new and x_removed

      for (int d = 0; d < D; ++d) {
        if (id_pblm[d]) {
          x_removed.push_back(x_unik[d]);
        } else {
          x_unik_new.push_back(x_unik[d]);
        }
      }
    } else {
      // x_unik_new is identical to x_unik
      x_unik_new = x_unik;
    }
  }
}
