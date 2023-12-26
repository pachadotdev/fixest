#include <cpp11.hpp>
#include <vector>

using namespace cpp11;
using namespace std;

// with cpp11 we write this
// considerations
// the original obsCluster is a "chunky" IntegerVector of dimension X*Y
// which is a matrix
[[cpp11::register]] list cpp11_get_fe_gnl(int Q, int N, doubles sumFE,
                                          integers_matrix<> dumMat,
                                          integers cluster_sizes,
                                          integers_matrix<> obsCluster) {
  int iter = 0, iterMax = 10000;
  int iter_loop = 0, iterMax_loop = 10000;
  int nb_coef = 0;

  writable::integers nb_ref(Q);  // nb_ref takes the nb of elements set as ref
  for (int q = 0; q < Q; q++) {
    // in rcpp: nb_coef += cluster_sizes(q);
    // in cpp11: we cannot use () for vectors, we use []
    nb_coef += cluster_sizes[q];
  }

  writable::doubles cluster_values(nb_coef);
  // vector<double> cluster_values(nb_coef);

  // IntegerVector cluster_visited(nb_coef);
  // integers cluster_visited(nb_coef);
  // this is not used! PROBLEM???

  std::vector<int *> pindex_cluster(Q);
  std::vector<int> index_cluster(nb_coef);
  for (int i = 0; i < nb_coef; i++) {
    index_cluster[i] = i;
  }
  pindex_cluster[0] = index_cluster.data();
  for (int q = 1; q < Q; q++) {
    // rcpp: pindex_cluster[q] = pindex_cluster[q - 1] + cluster_sizes(q - 1);
    // cpp11: we cannot use () for vectors, we use []
    pindex_cluster[q] = pindex_cluster[q - 1] + cluster_sizes[q - 1];
  }

  writable::integers start_cluster(nb_coef);
  writable::integers end_cluster(nb_coef);

  for (int i = 0; i < nb_coef; i++) {
    start_cluster[i] = 0;
    end_cluster[i] = 0;
  }

  int index;
  int k;

  for (int q = 0; q < Q; q++) {
    // rcpp: IntegerVector tableCluster(cluster_sizes(q));
    // cpp11: we cannot use () for vectors, we use []
    writable::integers tableCluster(cluster_sizes[q]);

    for (int i = 0; i < N; i++) {
      // cpp11: for a matrix we can extract entry Aij as A(i,j)
      k = dumMat(i, q);

      // rcpp: tableCluster(k) += 1;  // the number of obs per case
      // cpp11: we cannot use () for vectors, we use []
      tableCluster[k] += 1;
    }

    // rcpp: for (int k = 0; k < cluster_sizes(q); k++) {
    // cpp11: we cannot use () for vectors, we use []
    for (int k = 0; k < cluster_sizes[q]; k++) {
      int *pindex = pindex_cluster[q];
      index = pindex[k];

      // in rcpp this works, but in cpp11 this will create this error
      // error: ambiguous overload for ‘operator=’
      // (operand types are ‘cpp11::writable::r_vector<int>::proxy’ and
      // ‘cpp11::writable::r_vector<int>::proxy’)
      // start_cluster[index] = end_cluster[index - 1];
      // so we fill with 0s before hand and use += instead
      //   if (k == 0) {
      //     start_cluster[index] = 0;
      //     end_cluster[index] = tableCluster[k];
      //   } else {
      //     start_cluster[index] = end_cluster[index - 1];
      //     end_cluster[index] = end_cluster[index - 1] + tableCluster[k];
      //   }
      if (k == 0) {
        end_cluster[index] += tableCluster[k];
      } else {
        start_cluster[index] += end_cluster[index - 1];
        end_cluster[index] += end_cluster[index - 1] + tableCluster[k];
      }
    }
  }

  writable::integers_matrix<> mat_done(N, Q);

  // in rcpp we can fill with 0s with IntegerVector rowsums(N, 0)
  // in cpp11 we need to fill with a loop
  writable::integers rowsums(N);
  for (int i = 0; i < N; ++i) {
    rowsums[i] = 0;
  }

  writable::integers id2do(N);
  writable::integers id2do_next(N);

  // vector<int> id2do(N);
  // vector<int> id2do_next(N);

  int nb2do = N, nb2do_next = N;
  for (int i = 0; i < nb2do; i++) {
    id2do[i] = i;  // [] vs ()
    id2do_next[i] = i;
  }

  int qui_max, obs;
  int rs, rs_max;  // rs: row sum
  int id_cluster;
  double other_value;
  bool first;
  while (iter < iterMax) {
    iter++;
    if (iter == 1) {
      qui_max = 0;
    } else {
      qui_max = 0;
      rs_max = 0;
      for (int i = 0; i < nb2do; i++) {
        obs = id2do[i];
        rs = rowsums[obs];
        if (rs == Q - 2) {
          qui_max = obs;
          break;
        } else if (rs < Q && rs > rs_max) {
          qui_max = obs;
          rs_max = rs;
        } else if (qui_max == 0 && rs == 0) {
          qui_max = obs;
        }
      }
    }
    first = true;
    for (int q = 0; q < Q; q++) {
      // same consideration, in cpp11 Aij is A(i,j), not A[i,j]
      if (mat_done(qui_max, q) == 0) {
        if (first) {
          first = false;
        } else {
          id_cluster = dumMat(qui_max, q);
          int *pindex = pindex_cluster[q];
          index = pindex[id_cluster];
          cluster_values[index] = 0;
          for (int i = start_cluster[index]; i < end_cluster[index]; i++) {
            // error: no match for call to
            // ‘(cpp11::integers {aka cpp11::r_vector<int>}) (int&, int&)’
            // obs = obsCluster(i, q);
            // solution: change the class from "chunky vector" to matrix
            obs = obsCluster(i, q);
            mat_done(obs, q) = 1;
            rowsums[obs]++;
          }

          // error: no match for call to
          // ‘(cpp11::writable::integers {aka cpp11::writable::r_vector<int>})
          // (int&)’ nb_ref(q)++; solution: use [] instead of ()
          nb_ref[q]++;
        }
      }
    }
    iter_loop = 0;
    while (iter_loop < iterMax_loop) {
      iter_loop++;
      R_CheckUserInterrupt();
      if (iter_loop != 1) {
        nb2do = nb2do_next;
        for (int i = 0; i < nb2do; i++) {
          // WHY???????

          // error: ambiguous overload for ‘operator=’
          // (operand types are ‘cpp11::writable::r_vector<int>::proxy’ and
          // ‘cpp11::writable::r_vector<int>::proxy’)
          // id2do[i] = id2do_next[i];

          // this should resolve the ambiguity
          // id2do.at(i) = id2do_next.at(i);
          // not really, also error: ambiguous overload for ‘operator=’

          // maybe
          // id2do.operator[](i).operator=(id2do_next.operator[](i));
          // nope
          // error: call of overloaded
          // ‘operator=(cpp11::writable::r_vector<int>::proxy)’ is ambiguous

          // not either
          // int temp = id2do_next[i];
          // id2do[i] = temp;

          // ok, I gave up and used vector<int> instead of writable::integers
          // nope, it didn't work

          id2do[i] = 0 + id2do_next[i];
        }
      }
      nb2do_next = 0;
      for (int i = 0; i < nb2do; i++) {
        obs = id2do[i];
        rs = rowsums[obs];
        if (rs < Q - 1) {
          id2do_next[nb2do_next] = obs;
          nb2do_next++;
        } else if (rs == Q - 1) {
          int q;
          for (q = 0; mat_done(obs, q) != 0; q++) {
          }
          int *pindex = pindex_cluster[q];
          int index_select = pindex[dumMat(obs, q)];
          other_value = 0;
          for (int l = 0; l < Q; l++) {
            index = pindex_cluster[l][dumMat(obs, l)];
            // rcpp: other_value += cluster_values(index);
            // cpp11: we cannot use () for vectors, we use []
            other_value += cluster_values[index];
          }
          // rcpp: cluster_values(index_select) = sumFE(obs) - other_value;
          // cpp11: we cannot use () for vectors, we use []
          cluster_values[index_select] = sumFE[obs] - other_value;
          for (int j = start_cluster[index_select];
               j < end_cluster[index_select]; j++) {
            obs = obsCluster(j, q);
            mat_done(obs, q) = 1;
            rowsums[obs]++;
          }
        }
      }
      if (nb2do_next == nb2do) break;
    }
    if (iter_loop == iterMax_loop) {
      Rprintf(
          "Problem getting FE, maximum iterations reached (2nd order loop).");
    }
    if (nb2do_next == 0) break;
  }
  if (iter == iterMax) {
    Rprintf("Problem getting FE, maximum iterations reached (1st order loop).");
  }

  writable::list res(Q + 1);
  // maybe create a list and push back to it instead
  // writable::list res;

  int *pindex;
  for (int q = 0; q < Q; q++) {
    writable::doubles quoi(cluster_sizes[q]);
    // vector<double> quoi(cluster_sizes[q]);
    pindex = pindex_cluster[q];
    // rcpp: for (k = 0; k < cluster_sizes(q); k++) {
    // cpp11: we cannot use () for vectors, we use []
    for (k = 0; k < cluster_sizes[q]; k++) {
      index = pindex[k];
      // rcpp: quoi(k) = cluster_values(index);
      // error: ambiguous overload for ‘operator=’
      // (operand types are ‘cpp11::writable::r_vector<double>::proxy’ and
      // ‘cpp11::writable::r_vector<double>::proxy’)
      // quoi[k] = cluster_values[index];
      // not quoi.operator[](k).operator=(cluster_values.operator[](index));
      // either
      // error as well
      // double temp = cluster_values[index];
      // quoi[k] = temp;

      // ok, I gave up and used vector<double> instead of writable::doubles
      // nope, we need doubles to push to list
      quoi[k] = 0 + cluster_values[index];
    }

    // rcpp:: res(q) = quoi;
    // after using vector<int/double> instead of writable::integers/doubles
    // now it doesnt't like this
    res[q] = quoi;
    // error: no match for ‘operator=’ (operand types are ‘cpp11::list’ and
    // ‘std::vector<double>’)
    // push back maybe?
    // res.push_back({quoi});
  }
  // rcpp: res(Q) = nb_ref;
  //  this fails
  res[Q] = nb_ref;
  // push back maybe?
  // res.push_back({nb_ref});
  return res;
}
