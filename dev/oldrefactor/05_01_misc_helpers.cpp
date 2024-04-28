#include "05_0_misc.hpp"

[[cpp11::register]] std::string cpp_add_commas_(double x, int r = 1,
                                                bool whole = true) {
  // a bit like (but not exactly equal to) format(x, nsmall = 1, big.mark = ",")
  // but about 40-100 times faster for whole numbers => no trailing digits does
  // not accept vectors, although super easy to expand to vectors

  std::string x_str = std::to_string(static_cast<int>(abs(x)));
  std::string res;

  if (x < 0) {
    res.push_back('-');
    x = -x;
  }

  if (x < 1000) {
    res.insert(res.size(), x_str);
  } else {
    int n = x_str.size();
    int e = n; // e: exponent

    while (e > 0) {
      res.push_back(x_str[n - e]);
      --e;
      if (e > 1 && e % 3 == 0) {
        res.push_back(',');
      }
    }
  }

  double rest = x - floor(x);
  if ((rest != 0 || !whole) && r > 0) {
    // not a whole number

    res.push_back('.');

    if (r == 1) {
      res.push_back(std::to_string(static_cast<int>(round(rest * 10)))[0]);
    } else {
      double rounded_rest = round(rest * pow(10, r)) / pow(10, r);
      std::string rest_str = std::to_string(rounded_rest);
      int nr = rest_str.size();

      for (int i = 2; i < nr && i - 1 <= r; ++i) {
        res.push_back(rest_str[i]);
      }
    }
  }

  return res;
}

[[cpp11::register]] list cpp_find_never_always_treated_(integers cohort,
                                                        doubles period) {
  // Note that both cohort and period are sorted according to cohort

  writable::integers always_treated;
  // cohort_ref => the guys that will be out of the interaction
  writable::integers cohort_ref;

  int n = cohort.size();

  bool is_pos = false;
  bool is_neg = false;
  bool is_ok = false;

  int j = cohort[0];
  if (period[0] < 0) {
    is_neg = true;
  } else {
    is_pos = true;
  }

  for (int i = 1; i < n; ++i) {
    if (j == cohort[i]) {
      if (!is_ok) {
        // condition avoids checking period
        // only worth for very (very) long vectors

        if (period[i] < 0) {
          is_neg = true;
          is_ok = is_pos;
        } else {
          is_pos = true;
          is_ok = is_neg;
        }
      }
    } else {
      // we change IDs
      if (!is_ok) {
        if (is_pos)
          always_treated.push_back(j);
        cohort_ref.push_back(j);
      }

      // re-init
      j = cohort[i];
      is_neg = false;
      is_pos = false;
      is_ok = false;
    }
  }

  // Last element
  if (!is_ok) {
    if (is_pos)
      always_treated.push_back(j);
    cohort_ref.push_back(j);
  }

  writable::list res;
  res.push_back({"always_treated"_nm = always_treated});
  res.push_back({"ref"_nm = cohort_ref});

  return res;
}

[[cpp11::register]] integers cpp_get_first_item_(integers x, int n_items) {
  // observation id of the first occurrence
  // x ranges from 1 to n_items
  // we return indexes R style

  int n = x.size();

  writable::integers res(n_items);
  for (int i = 0; i < n_items; ++i) {
    res[i] = 0;
  }

  for (int i = 0; i < n; ++i) {
    if (res[x[i] - 1] == 0) {
      res[x[i] - 1] = i + 1;
    }
  }

  return res;
}

[[cpp11::register]] integers cpp_combine_clusters_(SEXP cluster_list,
                                                   integers index) {
  // cluster: list of integer vectors, each ranging from 1 to the number of
  // cases index: result of order() on the clusters

  if (TYPEOF(cluster_list) != VECSXP) {
    stop("Internal error: Only lists are accepted!");
  }

  int Q = Rf_length(cluster_list);

  int n = index.size();
  writable::integers res(n);

  // Loading the data
  vector<int *> pcluster(Q);
  for (int q = 0; q < Q; ++q) {
    SEXP cluster_q = VECTOR_ELT(cluster_list, q);

    pcluster[q] = INTEGER(cluster_q);
  }

  // the observation ID
  int obs = index[0] - 1;

  // vector holding the current value
  vector<int> current_value(Q);

  // initialization
  int counter = 1;
  res[obs] = counter;
  for (int q = 0; q < Q; ++q) {
    current_value[q] = pcluster[q][obs];
  }

  // we loop on the vector and flag values that are different
  int q = 0;
  for (int i = 1; i < n; ++i) {
    obs = index[i] - 1;

    for (q = 0; q < Q; ++q) {
      if (pcluster[q][obs] != current_value[q]) {
        break;
      }
    }

    // if the condition holds => means the values are different
    if (q < Q) {
      ++counter;
      // we save the new values
      for (; q < Q; ++q) {
        current_value[q] = pcluster[q][obs];
      }
    }

    res[obs] = counter;
  }

  return res;
}

[[cpp11::register]] list cpp_cut_(doubles x_sorted, doubles cut_points,
                                  integers is_included) {
  // x_sorted: no NA, sorted
  // cut_points: bounds
  // is_included: for each bound, if it is included or not
  //

  int N = x_sorted.size();
  int n_cuts = cut_points.size();

  bool is_int = true;
  for (int i = 0; i < N; ++i) {
    if (fabs(x_sorted[i] - round(x_sorted[i])) > 0.00000000001) {
      is_int = false;
      break;
    }
  }

  writable::integers x_int(N);
  for (int i = 0; i < N; i++) {
    x_int[i] = n_cuts + 1;
  }

  writable::integers isnt_empty(n_cuts + 1);
  writable::doubles value_min(n_cuts + 1);
  writable::doubles value_max(n_cuts + 1);

  int index = 0;
  bool first = true;
  bool include = is_included[0];
  double cutoff = cut_points[0];
  int i = 0;
  while (i < N) {
    if (include ? x_sorted[i] <= cutoff : x_sorted[i] < cutoff) {
      if (first) {
        isnt_empty[index] = true;
        value_min[index] = x_sorted[i];
        first = false;
      }

      x_int[i] = index + 1;
      ++i;
    } else {
      // we increment index
      // bins can be empty

      if (isnt_empty[index] && i > 0) {
        value_max[index] = x_sorted[i - 1];
      }

      ++index;

      if (index == n_cuts) {
        // last bin: we've done it at initialization
        isnt_empty[index] = true;
        value_min[index] = x_sorted[i];
        value_max[index] = x_sorted[N - 1];
        break;
      }

      include = is_included[index];
      cutoff = cut_points[index];
      first = true;
    }
  }

  if (index != n_cuts) {
    value_max[index] = x_sorted[N - 1];
  }

  writable::list res;
  res.push_back({"x_int"_nm = x_int});
  res.push_back({"isnt_empty"_nm = isnt_empty});
  res.push_back({"value_min"_nm = value_min});
  res.push_back({"value_max"_nm = value_max});
  res.push_back({"is_int"_nm = is_int});

  return res;
}

[[cpp11::register]] bool cpp_is_int_(SEXP x) {
  if (TYPEOF(x) == INTSXP) {
    return true;
  }

  if (TYPEOF(x) != REALSXP) {
    return false;
  }

  int N = Rf_length(x);
  double *px = REAL(x);

  bool is_int = true;
  for (int i = 0; i < N; ++i) {
    if (fabs(px[i] - round(px[i])) > 0.00000000001) {
      is_int = false;
      break;
    }
  }

  return is_int;
}

[[cpp11::register]] double cpp_hash_string_(std::string x) {
  // simple function hashing a string
  // used to identify tables in etable

  double res = std::hash<std::string>{}(x);

  return res;
}

[[cpp11::register]] bool cpp_isConstant_(doubles x) {
  // simple fun to see whether a variable is constant
  // it is unexpensive -- not the most useful function however:
  //		for 1e7 obs, you gain 50ms over var(x)==0, but it's still 50ms!

  int n = x.size();
  bool res = true;
  double value = x[0];
  for (int i = 1; i < n; i++) {
    if (x[i] != value) {
      res = false;
      break;
    }
  }

  return res;
}

[[cpp11::register]] bool cpp_any_na_null_(SEXP x) {
  // > twice faster than testing the two separately
  // x is a vector

  int n = Rf_length(x);
  double *px = REAL(x);

  for (int i = 0; i < n; ++i) {
    double x_tmp = px[i];
    if (std::isnan(x_tmp) || x_tmp == 0) {
      return true;
    }
  }

  return false;
}

[[cpp11::register]] list cpp_find_duplicates_(integers id, integers time) {
  // we check whether there are duplicated rows
  // if so, we provide information

  int n = id.size();
  int n_dup = 0;
  int obs_dup = 0;
  bool any_dup = false;

  /// finding the first duplicate value
  int i = 0;
  for (i = 1; i < n; ++i) {
    if (time[i - 1] == time[i]) {
      if (id[i - 1] == id[i]) {
        any_dup = true;
        break;
      }
    }
  }

  // if dup: we find out the nber
  if (any_dup) {
    obs_dup = i; // the 1 is implicitely added to make it r style
    int id_dup = id[i];
    int time_dup = time[i];
    n_dup = 2;
    while (++i < n && id_dup == id[i] && time_dup == time[i])
      n_dup++;
  }

  writable::list res;
  res.push_back({"n_dup"_nm = n_dup});
  res.push_back({"obs_dup"_nm = obs_dup});

  return res;
}

[[cpp11::register]] integers cpp_check_nested_(SEXP fe_list, SEXP cluster_list,
                                               integers fe_sizes, int n) {
  // Returns boolean vector of whether each FE is nested in the clusters

  int Q = Rf_length(fe_list);
  int G = Rf_length(cluster_list);

  // SEXP x0 = VECTOR_ELT(x, 0);

  writable::integers res(Q);

  for (int q = 0; q < Q; ++q) {
    int *pfe = INTEGER(VECTOR_ELT(fe_list, q));

    for (int g = 0; g < G; ++g) {
      vector<int> fe_clust(fe_sizes[q]);

      int *pclust = INTEGER(VECTOR_ELT(cluster_list, g));

      bool nested = true;
      int fe_value = 0;
      int clust_value = 0;
      for (int i = 0; i < n; ++i) {
        fe_value = pfe[i] - 1;
        clust_value = fe_clust[fe_value];
        if (clust_value == 0) {
          fe_clust[fe_value] = pclust[i];
        } else if (clust_value != pclust[i]) {
          nested = false;
          break;
        }
      }

      if (nested) {
        res[q] = 1;
        break;
      }
    }
  }

  return res;
}

[[cpp11::register]] int cpp_pgcd_(integers x) {
  // quick and dirty, but does not matter

  int n = x.size();

  if (n == 1) {
    return (x[0]);
  }

  bool ok = false;
  int pgcd = x[0];

  // the min
  for (int i = 1; i < n; ++i) {
    if (pgcd > x[i]) {
      pgcd = x[i];
    }
  }

  // the denom
  while (!ok && pgcd > 1) {
    ok = true;
    for (int i = 0; i < n; ++i) {
      if (x[i] % pgcd != 0) {
        pgcd--;
        ok = false;
        break;
      }
    }
  }

  return pgcd;
}

inline bool is_md_markup(const char *x, int i, int n) {
  if (x[i] != '*')
    return false;

  if (i + 1 >= n || x[i + 1] != '*')
    return true;
  if (i + 2 >= n || x[i + 2] != '*')
    return true;
  if (i + 3 < n && x[i + 3] == '*')
    return false;
  return true;
}

inline bool is_special_char(const char x) {
  return x == '&' || x == '%' || x == '_' || x == '^' || x == '#';
}

inline int find_id_markup(const char *x, int i, int n) {
  if (i + 1 >= n || x[i + 1] != '*')
    return 1;
  if (i + 2 >= n || x[i + 2] != '*')
    return 2;
  return 3;
}

std::string apply_escape_markup(const char *x) {
  // element 0 is the main character vector: it is never closed

  // int n = std::strlen(x); doesn't work with g++ but works with clang :S
  int n = strlen(x);
  std::vector<std::string> tmp_all(4, "");
  std::string tmp = "", res = "";

  std::vector<bool> is_open(4, false);
  std::vector<int> id_open;
  id_open.push_back(0);

  std::vector<std::string> open_tags{"", "\\textit{", "\\textbf{",
                                     "\\textbf{\\textit{"};
  std::vector<std::string> closing_tags{"", "}", "}", "}}"};
  std::vector<std::string> stars{"", "*", "**", "***"};

  int id_mkp = 0;

  int i = 0, i_save = 0;
  while (i < n) {
    if (x[i] == '$') {
      if (i > 0 && x[i - 1] == '\\') {
        tmp_all[id_mkp] += '$';
      } else {
        // This is an equation
        i_save = i;
        ++i;
        tmp = "$";
        while (i < n && (x[i] != '$' || x[i - 1] == '\\')) {
          tmp += x[i];
          ++i;
        }

        if (i == n) {
          // we went all the way without a closing $
          // => we come back, escape it and continue
          tmp_all[id_mkp] += "\\$";
          i = i_save;
        } else {
          tmp_all[id_mkp] += tmp + "$";
        }
      }
    } else if (is_special_char(x[i])) {
      if (i > 0 && x[i - 1] == '\\') {
        // we do nothing
        tmp_all[id_mkp] += x[i];
      } else {
        // we escape
        tmp_all[id_mkp] += "\\";
        tmp_all[id_mkp] += x[i];
      }
    } else if (x[i] == '\\') {
      if (i + 1 < n && x[i + 1] == '*') {
        // escape of MD markup
        ++i;
        while (i < n && x[i] == '*') {
          tmp_all[id_mkp] += '*';
          ++i;
        }
        --i; // compensated at the end of the main while
      } else {
        tmp_all[id_mkp] += '\\';
      }
    } else if (is_md_markup(x, i, n)) {
      id_mkp = find_id_markup(x, i, n);
      i += id_mkp - 1;

      if (is_open[id_mkp]) {
        for (int j = id_open.size() - 1; j >= 1; --j) {
          int id_j = id_open[j];
          int id_prev = id_open[j - 1];

          if (id_j == id_mkp) {
            // we apply the markup
            tmp = open_tags[id_mkp] + tmp_all[id_mkp] + closing_tags[id_mkp];

            tmp_all[id_prev] += tmp;

            is_open[id_mkp] = false;
            id_open.pop_back();
            tmp_all[id_mkp] = "";
            id_mkp = id_prev;
            break;
          } else {
            // we flush: the markup is invalid
            // example:  "** bonjour * les gens **"

            tmp = stars[id_j] + tmp_all[id_j];
            tmp_all[id_j] = "";
            tmp_all[id_prev] += tmp;

            is_open[id_j] = false;
            id_open.pop_back();
          }
        }
      } else {
        is_open[id_mkp] = true;
        id_open.push_back(id_mkp);
      }
    } else {
      tmp_all[id_mkp] += x[i];
    }
    ++i;
  }

  // we flush the rest if needed
  for (int j = id_open.size() - 1; j >= 1; --j) {
    int id_j = id_open[j];
    int id_prev = id_open[j - 1];

    tmp = stars[id_j] + tmp_all[id_j];
    tmp_all[id_j] = "";
    tmp_all[id_prev] += tmp;
  }

  return tmp_all[0];
}

[[cpp11::register]] strings cpp_escape_markup_(SEXP Rstr) {
  // cpp_escape_markup("**bonjour** *les* ***gens * \\***heureux*** ")
  // cpp_escape_markup("stars: 10%: *, 5%: **, 1%: ***")

  int n = LENGTH(Rstr);
  writable::strings res(n);

  if (n == 0) {
    return res;
  }

  for (int i = 0; i < n; ++i) {
    res[i] = apply_escape_markup(CHAR(STRING_ELT(Rstr, i)));
  }

  return res;
}
