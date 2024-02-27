#include "05_0_misc.hpp"

void initializePindexCluster(vector<int *> &pindex_cluster,
                             vector<int> &index_cluster,
                             const integers &cluster_sizes, int Q) {
  pindex_cluster[0] = index_cluster.data();
  for (int q = 1; q < Q; q++) {
    pindex_cluster[q] = pindex_cluster[q - 1] + cluster_sizes[q - 1];
  }
}

integers calculateTableCluster(const integers_matrix<> &dumMat,
                               const integers &cluster_sizes, int q, int N) {
  writable::integers tableCluster(cluster_sizes[q]);
  for (int i = 0; i < cluster_sizes[q]; i++) {
    tableCluster[i] = 0;
  }

  for (int i = 0; i < N; i++) {
    int k = dumMat(i, q);
    tableCluster[k] += 1;  // the number of obs per case
  }
  return tableCluster;
}

void setClusterIndices(int q, const integers &cluster_sizes,
                       const integers &tableCluster,
                       vector<int *> &pindex_cluster,
                       writable::integers &start_cluster,
                       writable::integers &end_cluster) {
  for (int k = 0; k < cluster_sizes[q]; k++) {
    int *pindex = pindex_cluster[q];
    int index = pindex[k];
    if (k == 0) {
      start_cluster[index] = 0;
      end_cluster[index] = tableCluster[k];
    } else {
      start_cluster[index] = static_cast<int>(end_cluster[index - 1]);
      end_cluster[index] = end_cluster[index - 1] + tableCluster[k];
    }
  }
}

int calculateQuiMax(const integers &id2do, const integers &rowsums, int nb2do,
                    int Q) {
  int qui_max = 0;
  int rs_max = -1;
  for (int i = 0; i < nb2do; i++) {
    int obs = id2do[i];
    int rs = rowsums[obs];
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
  return qui_max;
}

void updateId2do(int nb2do, writable::integers &id2do,
                 const integers &id2do_next) {
  // doesn't work with integers
  // copy(id2do_next.begin(), id2do_next.begin() + nb2do, id2do.begin());
  for (int i = 0; i < nb2do; ++i) {
    id2do[i] = id2do_next[i];
  }
}

void updateClusterValues(int obs, int Q, writable::integers_matrix<> &mat_done,
                         const integers_matrix<> &dumMat,
                         vector<int *> &pindex_cluster,
                         writable::doubles &cluster_values,
                         integers &start_cluster, integers &end_cluster,
                         integers_matrix<> &obsCluster,
                         writable::integers &rowsums, doubles &sumFE) {
  int q = 0;
  while (mat_done(obs, q) != 0) {
    q++;
  }
  int *pindex = pindex_cluster[q];
  int index_select = pindex[dumMat(obs, q)];
  double other_value = 0;
  for (int l = 0; l < Q; l++) {
    int index = pindex_cluster[l][dumMat(obs, l)];
    other_value += cluster_values[index];
  }
  cluster_values[index_select] = sumFE[obs] - other_value;
  for (int j = start_cluster[index_select]; j < end_cluster[index_select];
       j++) {
    int obs = obsCluster(j, q);
    mat_done(obs, q) = 1;
    rowsums[obs]++;
  }
}

[[cpp11::register]] list cpp_get_fe_gnl_(int Q, int N, doubles sumFE,
                                         integers_matrix<> dumMat,
                                         integers cluster_sizes,
                                         integers_matrix<> obsCluster) {
  int iter = 0, iterMax = 10000;
  int iter_loop = 0, iterMax_loop = 10000;
  int nb_coef = accumulate(cluster_sizes.begin(), cluster_sizes.end(), 0);

  writable::integers nb_ref(Q);
  for (int q = 0; q < Q; q++) {
    nb_ref[q] = 0;
  }

  writable::doubles cluster_values(nb_coef);
  for (int i = 0; i < nb_coef; i++) {
    cluster_values[i] = 0;
  }

  // writable::integers cluster_visited(nb_coef);
  // for (int i = 0; i < nb_coef; i++) {
  //   cluster_visited[i] = 0;
  // }

  vector<int *> pindex_cluster(Q);
  vector<int> index_cluster(nb_coef);
  iota(index_cluster.begin(), index_cluster.end(), 0);
  initializePindexCluster(pindex_cluster, index_cluster, cluster_sizes, Q);

  writable::integers start_cluster(nb_coef), end_cluster(nb_coef);
  for (int i = 0; i < nb_coef; i++) {
    start_cluster[i] = 0;
    end_cluster[i] = 0;
  }

  int index;
  for (int q = 0; q < Q; q++) {
    integers tableCluster = calculateTableCluster(dumMat, cluster_sizes, q, N);
    setClusterIndices(q, cluster_sizes, tableCluster, pindex_cluster,
                      start_cluster, end_cluster);
  }

  writable::integers_matrix<> mat_done(N, Q);
  for (int i = 0; i < N; i++) {
    for (int q = 0; q < Q; q++) {
      mat_done(i, q) = 0;
    }
  }

  writable::integers rowsums(N);
  for (int i = 0; i < N; i++) {
    rowsums[i] = 0;
  }

  writable::integers id2do(N);
  for (int i = 0; i < N; i++) {
    id2do[i] = i;
  }

  writable::integers id2do_next(N);
  for (int i = 0; i < N; i++) {
    id2do_next[i] = i;
  }

  int nb2do = N, nb2do_next = N;
  iota(id2do.begin(), id2do.end(), 0);
  iota(id2do_next.begin(), id2do_next.end(), 0);

  int qui_max, obs;
  int rs;  // rs: row sum
  int id_cluster;
  bool first;
  while (iter < iterMax) {
    iter++;
    if (iter == 1) {
      qui_max = 0;
    } else {
      qui_max = calculateQuiMax(id2do, rowsums, nb2do, Q);
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
          // VALGRIND ERROR: Conditional jump
          for (int i = start_cluster[index]; i < end_cluster[index]; i++) {
            obs = obsCluster(i, q);
            mat_done(obs, q) = 1;
            rowsums[obs]++;
          }
          nb_ref[q] += 1;
        }
      }
    }
    iter_loop = 0;
    while (iter_loop < iterMax_loop) {
      iter_loop++;
      check_user_interrupt();
      if (iter_loop != 1) {
        nb2do = nb2do_next;
        updateId2do(nb2do, id2do, id2do_next);
      }
      nb2do_next = 0;
      for (int i = 0; i < nb2do; i++) {
        obs = id2do[i];
        rs = rowsums[obs];
        if (rs < Q - 1) {
          id2do_next[nb2do_next] = obs;
          nb2do_next++;
        } else if (rs == Q - 1) {
          // VALGRIND ERROR: Conditional jump
          updateClusterValues(obs, Q, mat_done, dumMat, pindex_cluster,
                              cluster_values, start_cluster, end_cluster,
                              obsCluster, rowsums, sumFE);
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
  for (int q = 0; q < Q; q++) {
    auto &pindex = pindex_cluster[q];
    writable::doubles quoi(cluster_sizes[q]);
    for (int k = 0; k < cluster_sizes[q]; ++k) {
      quoi[k] = static_cast<double>(cluster_values[pindex[k]]);
    }
    res[q] = quoi;
  }
  res[Q] = nb_ref;

  return (res);
}
