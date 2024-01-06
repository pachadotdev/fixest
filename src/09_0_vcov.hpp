/*
__   __  ___   ___  __   __
\ \ / / / __| / _ \ \ \ / /
 \ V / | (__ | (_) | \ V /
  \_/   \___| \___/   \_/

Original Author: Laurent R. Berge
Refactored by Mauricio "Pacha" Vargas Sepulveda starting in Jun 2022
*/

#include "00_common.hpp"

#pragma once

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

// we define a class to access the data since access will be pretty intensive
class mat_row_scheme {
  int64_t K = 0;
  int64_t N = 0;
  int64_t n_total = 0;

 public:
  std::vector<double> mat;

  mat_row_scheme() = delete;
  mat_row_scheme(mat_row_scheme &);
  mat_row_scheme(doubles_matrix<> &);
  mat_row_scheme(int, int);

  double &operator()(int64_t i, int64_t k) { return mat[i * K + k]; }

  double *data() { return mat.data(); }

  void scale(double s) {
    for (int64_t i = 0; i < n_total; ++i) {
      mat[i] *= s;
    }
  }

  void check() {
    Rprintf("CHECK!\n");

    int cpt = 0;
    for (int64_t i = 0; i < N; ++i) {
      for (int64_t k = 0; k < K; ++k) {
        if (cpt++ <= 10)
          Rprintf("i = %d, k = %d, x = %d\n", i, k, mat[i * K + k]);
      }
    }
    Rprintf("...end\n");
  }

  int nrow() { return N; }
  int ncol() { return K; }
};

double to_sq(double x);
double dist_km(double lon_1, double lat_1, double cos_lat_1, double lon_2,
               double lat_2, double cos_lat_2);
double degree_to_radian(double x);
double fabs_lon(double x, double y);
double fabs_lat(double x, double y);
