#include <Rcpp.h>

#include <numeric>
#include <vector>

using namespace Rcpp;
using namespace std;

void initializePindexCluster(vector<int *> &pindex_cluster,
                             vector<int> &index_cluster,
                             const IntegerVector &cluster_sizes, int Q) {
  pindex_cluster[0] = index_cluster.data();
  for (int q = 1; q < Q; q++) {
    pindex_cluster[q] = pindex_cluster[q - 1] + cluster_sizes(q - 1);
  }
}

IntegerVector calculateTableCluster(const IntegerMatrix &dumMat,
                                    const IntegerVector &cluster_sizes, int q,
                                    int N) {
  IntegerVector tableCluster(cluster_sizes(q));
  for (int i = 0; i < N; i++) {
    int k = dumMat(i, q);
    tableCluster(k) += 1;  // the number of obs per case
  }
  return tableCluster;
}

void setClusterIndices(int q, const IntegerVector &cluster_sizes,
                       const IntegerVector &tableCluster,
                       vector<int *> &pindex_cluster,
                       IntegerVector &start_cluster,
                       IntegerVector &end_cluster) {
  for (int k = 0; k < cluster_sizes(q); k++) {
    int *pindex = pindex_cluster[q];
    int index = pindex[k];
    if (k == 0) {
      start_cluster[index] = 0;
      end_cluster[index] = tableCluster[k];
    } else {
      start_cluster[index] = end_cluster[index - 1];
      end_cluster[index] = end_cluster[index - 1] + tableCluster[k];
    }
  }
}

int calculateQuiMax(const Rcpp::IntegerVector &id2do,
                    const Rcpp::IntegerVector &rowsums, int nb2do, int Q) {
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

void updateId2do(int nb2do, Rcpp::IntegerVector &id2do,
                 const Rcpp::IntegerVector &id2do_next) {
  std::copy(id2do_next.begin(), id2do_next.begin() + nb2do, id2do.begin());
}

void updateClusterValues(int obs, int Q, IntegerMatrix &mat_done,
                         const IntegerMatrix &dumMat,
                         vector<int *> &pindex_cluster,
                         NumericVector &cluster_values,
                         IntegerVector &start_cluster,
                         IntegerVector &end_cluster, IntegerVector &obsCluster,
                         IntegerVector &rowsums, NumericVector &sumFE) {
  int q = 0;
  while (mat_done(obs, q) != 0) {
    q++;
  }
  int *pindex = pindex_cluster[q];
  int index_select = pindex[dumMat(obs, q)];
  double other_value = 0;
  for (int l = 0; l < Q; l++) {
    int index = pindex_cluster[l][dumMat(obs, l)];
    other_value += cluster_values(index);
  }
  cluster_values(index_select) = sumFE(obs) - other_value;
  for (int j = start_cluster[index_select]; j < end_cluster[index_select];
       j++) {
    int obs = obsCluster(j, q);
    mat_done(obs, q) = 1;
    rowsums[obs]++;
  }
}

// [[Rcpp::export]]
List rcpp_get_fe_gnl(int Q, int N, NumericVector sumFE, IntegerMatrix dumMat,
                     IntegerVector cluster_sizes, IntegerVector obsCluster) {
  int iter = 0, iterMax = 10000;
  int iter_loop = 0, iterMax_loop = 10000;
  IntegerVector nb_ref(Q);  // nb_ref takes the nb of elements set as ref
  int nb_coef = accumulate(cluster_sizes.begin(), cluster_sizes.end(), 0);
  NumericVector cluster_values(nb_coef);
  IntegerVector cluster_visited(nb_coef);
  vector<int *> pindex_cluster(Q);
  vector<int> index_cluster(nb_coef);
  iota(index_cluster.begin(), index_cluster.end(), 0);
  initializePindexCluster(pindex_cluster, index_cluster, cluster_sizes, Q);
  IntegerVector start_cluster(nb_coef), end_cluster(nb_coef);
  int index;
  int k;
  for (int q = 0; q < Q; q++) {
    IntegerVector tableCluster =
        calculateTableCluster(dumMat, cluster_sizes, q, N);
    setClusterIndices(q, cluster_sizes, tableCluster, pindex_cluster,
                      start_cluster, end_cluster);
  }
  IntegerMatrix mat_done(N, Q);
  IntegerVector rowsums(N, 0);
  IntegerVector id2do(N);
  IntegerVector id2do_next(N);
  int nb2do = N, nb2do_next = N;
  iota(id2do.begin(), id2do.end(), 0);
  iota(id2do_next.begin(), id2do_next.end(), 0);
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
      int qui_max = calculateQuiMax(id2do, rowsums, nb2do, Q);
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
          cluster_values(index) = 0;
          for (int i = start_cluster[index]; i < end_cluster[index]; i++) {
            obs = obsCluster(i, q);
            mat_done(obs, q) = 1;
            rowsums[obs]++;
          }
          nb_ref(q)++;
        }
      }
    }
    iter_loop = 0;
    while (iter_loop < iterMax_loop) {
      iter_loop++;
      R_CheckUserInterrupt();
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
  List res(Q + 1);
  for (int q = 0; q < Q; q++) {
    auto &pindex = pindex_cluster[q];
    NumericVector quoi(cluster_sizes[q]);
    for (int k = 0; k < cluster_sizes[q]; ++k) {
      quoi[k] = cluster_values[pindex[k]];
    }
    res[q] = quoi;
  }
  res(Q) = nb_ref;
  return (res);
}
