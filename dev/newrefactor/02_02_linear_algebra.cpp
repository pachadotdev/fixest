// We introduce a class that handles varying types of SEXP and behaves as a
// regular matrix

// Later I may use the number of observations
// I don't at the moment because everything is strongly checked beforehand

#include "02_0_demeaning.hpp"

sVec::sVec(SEXP x) {
  if (TYPEOF(x) == REALSXP) {
    is_int = false;
    p_dble = REAL(x);
  } else if (TYPEOF(x) == INTSXP) {
    is_int = true;
    p_int = INTEGER(x);
  } else {
    stop("The current SEXP type is not supported by the sVec class.");
  }
}

sMat::sMat(SEXP x, bool single_obs = false) {
  if (TYPEOF(x) == VECSXP) {
    // x can be a list of either vectors or matrices

    int L = Rf_length(x);

    for (int l = 0; l < L; ++l) {
      SEXP xx = VECTOR_ELT(x, l);
      SEXP dim = Rf_getAttrib(xx, R_DimSymbol);

      int n_tmp = 0, K_tmp = 0;

      if (Rf_length(dim) == 0) {
        // vector
        n_tmp = Rf_length(xx);
        K_tmp = 1;
      } else {
        int *pdim = INTEGER(dim);
        n_tmp = pdim[0];
        K_tmp = pdim[1];
      }

      // we set the number of rows at the first iteration
      if (l == 0) {
        n = n_tmp;
      } else {
        if (n != n_tmp)
          stop(
              "When setting up the class sMat: The number of observations in "
              "the list is not coherent across columns.");
      }

      K += K_tmp;

      if (TYPEOF(xx) == REALSXP) {
        double *p_x = REAL(xx);
        for (int k = 0; k < K_tmp; ++k) {
          p_sVec.push_back(sVec(p_x));
          if (k + 1 < K_tmp) p_x += n;
        }

      } else if (TYPEOF(xx) == INTSXP) {
        int *p_x = INTEGER(xx);
        for (int k = 0; k < K_tmp; ++k) {
          p_sVec.push_back(sVec(p_x));
          if (k + 1 < K_tmp) p_x += n;
        }
      } else {
        stop("The current SEXP type is not supported by the sMat class.");
      }
    }

  } else {
    // Matrix or vector

    SEXP dim = Rf_getAttrib(x, R_DimSymbol);

    if (Rf_length(dim) == 0) {
      // vector
      n = Rf_length(x);
      K = 1;
    } else {
      const int *pdim = INTEGER(dim);
      n = pdim[0];
      K = pdim[1];
    }

    if (!single_obs && (n == 1 && K == 1)) {
      // => absence of data
      n = 0;
      K = 0;

    } else if (TYPEOF(x) == REALSXP) {
      double *p_x = REAL(x);
      for (int k = 0; k < K; ++k) {
        p_sVec.push_back(sVec(p_x));
        if (k + 1 < K) p_x += n;
      }

    } else if (TYPEOF(x) == INTSXP) {
      int *p_x = INTEGER(x);
      for (int k = 0; k < K; ++k) {
        p_sVec.push_back(sVec(p_x));
        if (k + 1 < K) p_x += n;
      }
    } else {
      stop("The current SEXP type is not supported by the sMat class.");
    }
  }
}

sVec sMat::operator[](int k) { return p_sVec[k]; }

double sMat::operator()(int i, int k) { return p_sVec[k][i]; }

inline double &simple_mat_with_id::operator()(int id, int i, int j) {
  if (id != id_current) {
    id_current = id;
    px_current = px0 + n_total * id;
  }

  return px_current[i + nrow * j];
}

inline double &simple_mat_with_id::operator()(int id, int i) {
  if (id != id_current) {
    id_current = id;
    px_current = px0 + n_total * id;
  }

  return px_current[i];
}
