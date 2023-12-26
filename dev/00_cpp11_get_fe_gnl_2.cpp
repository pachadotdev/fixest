#include <cpp11.hpp>
#include <vector>

using namespace cpp11;
using namespace std;

// with cpp11 we write this
// considerations
// the original obsCluster is a "chunky" IntegerVector of dimension X*Y
// which is a matrix
[[cpp11::register]] list cpp11_get_fe_gnl(int Q, int N, writable::doubles sumFE,
                                          integers_matrix<> dumMat,
                                          integers cluster_sizes,
                                          integers_matrix<> obsCluster) {
  int iter = 0, iterMax = 10000;
  int iter_loop = 0, iterMax_loop = 10000;
  int nb_coef = 0;
  writable::integers nb_ref(Q);  // nb_ref takes the nb of elements set as ref
  for (int q = 0; q < Q; q++) {
    nb_coef += cluster_sizes[q];
  }
  writable::doubles cluster_values(nb_coef);
  writable::integers cluster_visited(nb_coef);
  std::vector<int *> pindex_cluster(Q);
  std::vector<int> index_cluster(nb_coef);
  for (int i = 0; i < nb_coef; i++) {
    index_cluster[i] = i;
  }
  pindex_cluster[0] = index_cluster.data();
  for (int q = 1; q < Q; q++) {
    pindex_cluster[q] = pindex_cluster[q - 1] + cluster_sizes[q - 1];
  }

  writable::integers start_cluster(nb_coef), end_cluster(nb_coef);
  for (int i = 0; i < nb_coef; ++i) {
    start_cluster[i] = 0;
    end_cluster[i] = 0;
  }

  int index;
  int k;
  for (int q = 0; q < Q; q++) {
    writable::integers tableCluster(cluster_sizes[q]);
    for (int i = 0; i < N; i++) {
      k = dumMat(i, q);
      tableCluster[k] += 1;  // the number of obs per case
    }
    for (int k = 0; k < cluster_sizes[q]; k++) {
      int *pindex = pindex_cluster[q];
      index = pindex[k];
      // index = start(q) + k;
      if (k == 0) {
        start_cluster[index] = 0;
        // end_cluster[index] = tableCluster[k];
        end_cluster[index] = static_cast<int>(tableCluster[k]);
      } else {
        // start_cluster[index] = end_cluster[index - 1];
        start_cluster[index] = static_cast<int>(end_cluster[index - 1]);

        end_cluster[index] = end_cluster[index - 1] + tableCluster[k];
      }
    }
  }

  writable::integers_matrix<> mat_done(N, Q);
  // fill matrix with 0s
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < Q; ++j) {
      mat_done(i, j) = 0;
    }
  }

  // writable::integers rowsums(N, 0);
  // in rcpp we can fill with 0s with IntegerVector rowsums(N, 0)
  // in cpp11 we need to fill with a loop
  writable::integers rowsums(N);
  for (int i = 0; i < N; ++i) {
    rowsums[i] = 0;
  }

  writable::integers id2do(N);
  writable::integers id2do_next(N);
  int nb2do = N, nb2do_next = N;
  for (int i = 0; i < nb2do; i++) {
    id2do[i] = i;
    id2do_next[i] = i;
  }
  int qui_max = 0, obs;
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
      if (mat_done(qui_max, q) == 0) {
        if (first) {
          first = false;
        } else {
          id_cluster = dumMat(qui_max, q);
          int *pindex = pindex_cluster[q];
          index = pindex[id_cluster];
          cluster_values[index] = 0;
          for (int i = start_cluster[index]; i < end_cluster[index]; i++) {
            obs = obsCluster(i, q);
            mat_done(obs, q) = 1;
            rowsums[obs]++;
          }
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
          // id2do[i] = id2do_next[i];
          id2do[i] = static_cast<int>(id2do_next[i]);
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
            other_value += cluster_values[index];
          }
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
      Rprintf("Problem getting FE");
    }
    if (nb2do_next == 0) break;
  }
  if (iter == iterMax) {
    Rprintf("Problem getting FE");
  }
  writable::list res(Q + 1);
  int *pindex;
  for (int q = 0; q < Q; q++) {
    writable::doubles quoi(cluster_sizes[q]);
    pindex = pindex_cluster[q];
    for (k = 0; k < cluster_sizes[q]; k++) {
      index = pindex[k];
      // quoi[k] = cluster_values[index];
      quoi[k] = static_cast<double>(cluster_values[index]);
    }
    res[q] = quoi;
  }
  res[Q] = nb_ref;
  return (res);
}