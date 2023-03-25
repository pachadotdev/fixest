#include "04_0_linear_model.h"

// to print elements with cout << X(i,k); we use
// #include <iostream>
// using namespace std;

void mp_sparse_XtX(writable::doubles_matrix<> &XtX, const vector<int> &n_j, const vector<int> &start_j, const vector<int> &all_i, const vector<double> &x, const doubles_matrix<> &X, int nthreads)
{

    int K = X.ncol();

    #pragma omp parallel for num_threads(nthreads) schedule(static, 1)
    for (int j1 = 0; j1 < K; ++j1)
    {

        int k = 0, start = 0, end = 0;

        // Upper triangle
        for (int j2 = j1; j2 < K; ++j2)
        {

            if (n_j[j1] < n_j[j2])
            {
                start = start_j[j1];
                end = start_j[j1 + 1];
                k = j2;
            }
            else
            {
                start = start_j[j2];
                end = start_j[j2 + 1];
                k = j1;
            }

            double value = 0;
            for (int index = start; index < end; ++index)
            {
                value += X(all_i[index], k) * x[index];
            }

            if (value == 0)
                continue;

            XtX(j1, j2) = value;
            XtX(j2, j1) = value;
        }
    }
}

void mp_sparse_Xty(writable::doubles &Xty, const vector<int> &start_j, const vector<int> &all_i, const vector<double> &x, const double *y, int nthreads)
{

    int K = Xty.size();

#pragma omp parallel for num_threads(nthreads)
    for (int j = 0; j < K; ++j)
    {

        int start = start_j[j];
        int end = start_j[j + 1];

        double value = 0;
        for (int index = start; index < end; ++index)
        {
            value += y[all_i[index]] * x[index];
        }

        if (value == 0)
            continue;

        Xty[j] = value;
    }
}

void mp_XtX(writable::doubles_matrix<> &XtX, const doubles_matrix<> &X, const doubles_matrix<> &wX, int nthreads)
{

    int N = X.nrow();
    int K = X.ncol();

    // specific scheme for large N and K == 1
    if (K == 1)
    {

        vector<double> all_values(nthreads, 0);
        vector<int> bounds = set_parallel_scheme(N, nthreads);

#pragma omp parallel for num_threads(nthreads)
        for (int t = 0; t < nthreads; ++t)
        {
            double val = 0;
            for (int i = bounds[t]; i < bounds[t + 1]; ++i)
            {
                val += X(i, 0) * wX(i, 0);
            }
            all_values[t] = val;
        }

        double value = 0;
        for (int t = 0; t < nthreads; ++t)
        {
            value += all_values[t];
        }

        XtX(0, 0) = value;
    }
    else
    {
        // We use this trick to even out the load on the threads
        int nValues = K * (K + 1) / 2;
        vector<int> all_i, all_j;
        for (int i = 0; i < K; ++i)
        {
            for (int j = i; j < K; ++j)
            {
                all_i.push_back(i);
                all_j.push_back(j);
            }
        }

#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
        for (int index = 0; index < nValues; ++index)
        {
            int k_row = all_i[index];
            int k_col = all_j[index];

            double val = 0;
            for (int i = 0; i < N; ++i)
            {
                val += X(i, k_row) * wX(i, k_col);
            }

            XtX(k_row, k_col) = val;
            XtX(k_col, k_row) = val;
        }
    }
}

void mp_Xty(writable::doubles &Xty, const doubles_matrix<> &X, const double *y, int nthreads)
{

    int N = X.nrow();
    int K = X.ncol();

    if (K == 1)
    {

        vector<double> all_values(nthreads, 0);
        vector<int> bounds = set_parallel_scheme(N, nthreads);

        #pragma omp parallel for num_threads(nthreads)
        for (int t = 0; t < nthreads; ++t)
        {
            double val = 0;
            for (int i = bounds[t]; i < bounds[t + 1]; ++i)
            {
                val += X(i, 0) * y[i];
            }
            all_values[t] = val;
        }

        double value = 0;
        for (int t = 0; t < nthreads; ++t)
        {
            value += all_values[t];
        }

        Xty[0] = value;
    }
    else
    {

    #pragma omp parallel for num_threads(nthreads)
        for (int j = 0; j < K; ++j)
        {
            double val = 0;
            for (int i = 0; i < N; ++i)
            {
                val += X(i, j) * y[i];
            }

            Xty[j] = val;
        }
    }
}

[[cpp11::register]] list cpp_sparse_products_(doubles_matrix<> X, doubles w, SEXP y, bool correct_0w = false, int nthreads = 1)
{

    int N = X.nrow();
    int K = X.ncol();

    bool isWeight = w.size() > 1;

    bool is_y_list = TYPEOF(y) == VECSXP;

    writable::doubles_matrix<> XtX(K, K);

    if (sparse_check(X) == false)
    {
        // NOT SPARSE

        writable::list res;

        // TODO: probably there's a much better way to obtain wX

        writable::doubles_matrix<> wX(N,K);

        if (isWeight)
        {
            #pragma omp parallel for num_threads(nthreads)
            for (int k = 0; k < K; ++k)
            {
                for (int i = 0; i < N; ++i)
                {
                    wX(i, k) = w[i] * X(i, k);
                }
            }
        } else {
            #pragma omp parallel for num_threads(nthreads)
            for (int k = 0; k < K; ++k)
            {
                for (int i = 0; i < N; ++i)
                {
                    wX(i, k) = X(i, k);
                }
            }
        }

        // XtX
        mp_XtX(XtX, X, wX, nthreads);
        res.push_back({"XtX"_nm = XtX});

        // Xty
        if (is_y_list)
        {
            int n_vars_y = Rf_length(y);
            writable::list Xty(n_vars_y);

            for (int v = 0; v < n_vars_y; ++v)
            {
                writable::doubles Xty_tmp(K);
                mp_Xty(Xty_tmp, wX, REAL(VECTOR_ELT(y, v)), nthreads);
                Xty[v] = Xty_tmp;
            }

            res.push_back({"Xty"_nm = Xty});
        }
        else
        {
            writable::doubles Xty(K);
            mp_Xty(Xty, wX, REAL(y), nthreads);
            res.push_back({"Xty"_nm = Xty});
        }

        return res;
    }

    //
    // SPARSE case
    //

    vector<int> n_j(K, 0);
    vector<int> start_j(K + 1, 0);
    vector<int> all_i;
    vector<double> x;

    set_sparse(n_j, start_j, all_i, x, X, w);

    writable::list res;

    // XtX
    mp_sparse_XtX(XtX, n_j, start_j, all_i, x, X, nthreads);
    res.push_back({"XtX"_nm = XtX});

    // Xty
    if (is_y_list)
    {
        int n_vars_y = Rf_length(y);
        writable::list Xty(n_vars_y);

        for (int v = 0; v < n_vars_y; ++v)
        {
            writable::doubles Xty_tmp(K);
            mp_sparse_Xty(Xty_tmp, start_j, all_i, x, REAL(VECTOR_ELT(y, v)), nthreads);
            Xty[v] = Xty_tmp;
        }

        res.push_back({"Xty"_nm = Xty});
    }
    else
    {
        writable::doubles Xty(K);
        mp_sparse_Xty(Xty, start_j, all_i, x, REAL(y), nthreads);
        res.push_back({"Xty"_nm = Xty});
    }

    return res;
}

[[cpp11::register]] doubles_matrix<> cpppar_crossprod_(doubles_matrix<> X, doubles w, int nthreads)
{

    int N = X.nrow();
    int K = X.ncol();

    bool isWeight = false;
    if (w.size() > 1)
    {
        isWeight = true;
    }

    writable::doubles_matrix<> XtX(K, K);

    if (sparse_check(X) == false)
    {
        // NOT SPARSE

        if (isWeight)
        {
            // I don't use the sqrt, because I use the function when weights are negative too (ll_d2 used as 'weight')
            writable::doubles_matrix<> wX(N,K);
            for (int k = 0; k < K; ++k)
            {
                for (int i = 0; i < N; ++i)
                {
                    wX(i, k) = X(i, k) * w[i];
                }
            }

            // XtX
            mp_XtX(XtX, X, wX, nthreads);
        }
        else
        {
            // Identical, but no need to make a copy of X or y which can be expensive

            // XtX
            mp_XtX(XtX, X, X, nthreads);
        }

        return XtX;
    }

    //
    // SPARSE case
    //

    vector<int> n_j(K, 0);
    vector<int> start_j(K + 1, 0);
    vector<int> all_i;
    vector<double> x;

    set_sparse(n_j, start_j, all_i, x, X, w);

    // XtX
    mp_sparse_XtX(XtX, n_j, start_j, all_i, x, X, nthreads);

    return XtX;
}

[[cpp11::register]] doubles_matrix<> cpp_mat_reconstruct_(writable::doubles_matrix<> X, logicals id_excl)
{

    int K = id_excl.size();
    int K_small = X.ncol();

    writable::doubles_matrix<> res(K, K);

    int n_col_excl = 0;
    for (int j = 0; j < K_small; ++j)
    {

        while (id_excl[j + n_col_excl])
            ++n_col_excl;

        int col = j + n_col_excl;
        int n_row_excl = 0;
        for (int i = 0; i < K_small; ++i)
        {
            while (id_excl[i + n_row_excl])
                ++n_row_excl;
            res(i + n_row_excl, col) += X(i, j);
        }
    }

    return res;
}

void mp_ZXtZX(writable::doubles_matrix<> &ZXtZX, const doubles_matrix<> &XtX, const doubles_matrix<> &X, const doubles_matrix<> &Z, const doubles_matrix<> &wZ, int nthreads)
{

    int N = Z.nrow();
    int K1 = Z.ncol();

    bool isX = X.nrow() > 1;
    int K2 = isX ? X.ncol() : 0;

    // First: we copy XtX
    for (int k = 0; k < K2; ++k)
    {
        for (int l = 0; l < K2; ++l)
        {
            ZXtZX(k + K1, l + K1) = XtX(k, l);
        }
    }

    //
    // The K2 x K1 band
    //

    // l: index of Z
    // k: index of X

    int nValues = K2 * K1;
    vector<int> all_l, all_k;
    for (int l = 0; l < K1; ++l)
    {
        for (int k = 0; k < K2; ++k)
        {
            all_l.push_back(l);
            all_k.push_back(k);
        }
    }

    #pragma omp parallel for num_threads(nthreads)
    for (int index = 0; index < nValues; ++index)
    {
        int l = all_l[index];
        int k = all_k[index];

        double val = 0;
        for (int i = 0; i < N; ++i)
        {
            val += X(i, k) * wZ(i, l);
        }

        ZXtZX(K1 + k, l) = val;
        ZXtZX(l, K1 + k) = val;
    }

    //
    // The K1 * K1 mat
    //

    // k, l: indexes of Z
    nValues = K1 * (K1 + 1) / 2;
    all_l.clear();
    all_k.clear();
    for (int l = 0; l < K1; ++l)
    {
        for (int k = l; k < K1; ++k)
        {
            all_l.push_back(l);
            all_k.push_back(k);
        }
    }

    #pragma omp parallel for num_threads(nthreads)
    for (int index = 0; index < nValues; ++index)
    {
        int l = all_l[index];
        int k = all_k[index];

        double val = 0;
        for (int i = 0; i < N; ++i)
        {
            val += Z(i, k) * wZ(i, l);
        }

        ZXtZX(k, l) = val;
        ZXtZX(l, k) = val;
    }
}

void mp_ZXtu(writable::doubles &ZXtu, const doubles_matrix<> &X, const doubles_matrix<> &Z, const double *u, int nthreads)
{

    int N = Z.nrow();
    int K1 = Z.ncol();

    bool isX = X.nrow() > 1;
    int K2 = isX ? X.ncol() : 0;

    #pragma omp parallel for num_threads(nthreads)
    for (int k = 0; k < (K2 + K1); ++k)
    {
        double val = 0;
        for (int i = 0; i < N; ++i)
        {
            if (k < K1)
            {
                val += Z(i, k) * u[i];
            }
            else
            {
                val += X(i, k - K1) * u[i];
            }
        }

        ZXtu[k] = val;
    }
}

void mp_sparse_ZXtZX(writable::doubles_matrix<> &ZXtZX, const doubles_matrix<> &XtX, const vector<int> &n_j, const vector<int> &start_j, const vector<int> &all_i, const vector<double> &x, const doubles_matrix<> &X, const doubles_matrix<> &Z, const doubles_matrix<> &wZ, int nthreads)
{

    int N = Z.nrow();
    int K1 = Z.ncol();

    bool isX = X.nrow() > 1;
    int K2 = isX ? X.ncol() : 0;

    // First: we copy XtX
    for (int k = 0; k < K2; ++k)
    {
        for (int l = 0; l < K2; ++l)
        {
            ZXtZX(k + K1, l + K1) = XtX(k, l);
        }
    }

    // the K2 x K1 band
    for (int l = 0; l < K1; ++l)
    {
        #pragma omp parallel for num_threads(nthreads) schedule(static, 1)
        for (int j = 0; j < K2; ++j)
        {

            int start = start_j[j];
            int end = start_j[j + 1];

            double value = 0;
            for (int index = start; index < end; ++index)
            {
                value += Z(all_i[index], l) * x[index];
            }

            ZXtZX(j + K1, l) = value;
            ZXtZX(l, j + K1) = value;
        }
    }

    //
    // The K1 * K1 mat
    //

    // k, l: indexes of Z
    int nValues = K1 * (K1 + 1) / 2;
    vector<int> all_l, all_k;
    for (int l = 0; l < K1; ++l)
    {
        for (int k = l; k < K1; ++k)
        {
            all_l.push_back(l);
            all_k.push_back(k);
        }
    }

    #pragma omp parallel for num_threads(nthreads)
    for (int index = 0; index < nValues; ++index)
    {
        int l = all_l[index];
        int k = all_k[index];

        double val = 0;
        for (int i = 0; i < N; ++i)
        {
            val += Z(i, k) * wZ(i, l);
        }

        ZXtZX(k, l) = val;
        ZXtZX(l, k) = val;
    }
}

void mp_sparse_ZXtu(writable::doubles &ZXtu, const vector<int> &start_j, const vector<int> &all_i, const vector<double> &x, const double *u, const doubles_matrix<> &X, const doubles_matrix<> &wZ, int nthreads)
{

    int N = wZ.nrow();
    int K1 = wZ.ncol();

    bool isX = X.nrow() > 1;
    int K2 = isX ? X.ncol() : 0;

    #pragma omp parallel for num_threads(nthreads)
    for (int k = 0; k < (K2 + K1); ++k)
    {

        double value = 0;
        if (k < K1)
        {
            for (int i = 0; i < N; ++i)
            {
                value += u[i] * wZ(i, k);
            }
        }
        else
        {
            int start = start_j[k - K1];
            int end = start_j[k - K1 + 1];

            for (int index = start; index < end; ++index)
            {
                value += u[all_i[index]] * x[index];
            }
        }

        ZXtu[k] = value;
    }
}
