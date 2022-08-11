#include "04_linear_model.h"

[[cpp11::register]] list cpp_iv_products(doubles_matrix<> X, SEXP y, doubles_matrix<> Z, SEXP u, doubles w, int nthreads)
{
    // We compute the following:
    // - X'X
    // - X'y
    // - (ZX)'(ZX)
    // - (ZX)'u
    //
    // Z: IV mat, u: endo regs

    // Note that X can be equal to 0 (ie no X)

    int N = Z.nrow();
    int K1 = Z.ncol();

    bool isX = X.nrow() > 1;
    int K2 = isX ? X.ncol() : 0;

    bool isWeight = w.size() > 1;

    bool is_y_list = TYPEOF(y) == VECSXP;

    writable::doubles_matrix<> XtX(K2, K2);
    writable::doubles_matrix<> ZXtZX(K2 + K1, K2 + K1);

    writable::doubles_matrix<> wZ(Z);
    if (isWeight)
    {
        for (int k = 0; k < K1; ++k)
        {
            for (int i = 0; i < N; ++i)
            {
                wZ(i, k) *= w[i];
            }
        }
    }

    if (sparse_check(X) == false)
    {
        // NOT SPARSE

        writable::list res;
        writable::doubles_matrix<> wX(X);
        if (isWeight)
        {
            for (int k = 0; k < K2; ++k)
            {
                for (int i = 0; i < N; ++i)
                {
                    wX(i, k) *= w[i];
                }
            }
        }

        // XtX
        mp_XtX(XtX, X, wX, nthreads);
        res["XtX"] = XtX;

        // (ZX)'(ZX)
        mp_ZXtZX(ZXtZX, XtX, X, Z, wZ, nthreads);
        res["ZXtZX"] = ZXtZX;

        // Xty
        if (!isX)
        {
            writable::doubles Xty(1);
            res["Xty"] = Xty;
        }
        else if (is_y_list)
        {
            int n_vars_y = Rf_length(y);
            writable::list Xty(n_vars_y);

            for (int v = 0; v < n_vars_y; ++v)
            {
                writable::doubles Xty_tmp(K2);
                mp_Xty(Xty_tmp, wX, REAL(VECTOR_ELT(y, v)), nthreads);
                Xty[v] = Xty_tmp;
            }

            res["Xty"] = Xty;
        }
        else
        {
            writable::doubles Xty(K2);
            mp_Xty(Xty, wX, REAL(y), nthreads);
            res["Xty"] = Xty;
        }

        // (ZX)'u
        // u is always a list
        int n_vars_u = Rf_length(u);
        writable::list ZXtu(n_vars_u);

        for (int v = 0; v < n_vars_u; ++v)
        {
            writable::doubles ZXtu_tmp(K2 + K1);
            mp_ZXtu(ZXtu_tmp, wX, wZ, REAL(VECTOR_ELT(u, v)), nthreads);
            ZXtu[v] = ZXtu_tmp;
        }

        res["ZXtu"] = ZXtu;

        return res;
    }

    //
    // SPARSE case
    //

    vector<int> n_j(K2 + !isX, 0);
    vector<int> start_j(K2 + !isX + 1, 0);
    vector<int> all_i;
    vector<double> x;

    set_sparse(n_j, start_j, all_i, x, X, w);

    writable::list res;

    // XtX
    mp_sparse_XtX(XtX, n_j, start_j, all_i, x, X, nthreads);
    res["XtX"] = XtX;

    // (xZ)'(ZX)
    mp_sparse_ZXtZX(ZXtZX, XtX, n_j, start_j, all_i, x, X, Z, wZ, nthreads);
    res["ZXtZX"] = ZXtZX;

    // Xty
    if (!isX)
    {
        writable::doubles Xty(1);
        res["Xty"] = Xty;
    }
    else if (is_y_list)
    {
        int n_vars_y = Rf_length(y);
        writable::list Xty(n_vars_y);

        for (int v = 0; v < n_vars_y; ++v)
        {
            writable::doubles Xty_tmp(K2);
            mp_sparse_Xty(Xty_tmp, start_j, all_i, x, REAL(VECTOR_ELT(y, v)), nthreads);
            Xty[v] = Xty_tmp;
        }

        res["Xty"] = Xty;
    }
    else
    {
        writable::doubles Xty(K2);
        mp_sparse_Xty(Xty, start_j, all_i, x, REAL(y), nthreads);
        res["Xty"] = Xty;
    }

    // (ZX)'u
    int n_vars_u = Rf_length(u);
    writable::list ZXtu(n_vars_u);

    for (int v = 0; v < n_vars_u; ++v)
    {
        writable::doubles ZXtu_tmp(K2 + K1);
        mp_sparse_ZXtu(ZXtu_tmp, start_j, all_i, x, REAL(VECTOR_ELT(u, v)), X, wZ, nthreads);
        ZXtu[v] = ZXtu_tmp;
    }

    res["ZXtu"] = ZXtu;

    return res;
}

[[cpp11::register]] list cpp_iv_product_completion(doubles_matrix<> XtX, doubles Xty, doubles_matrix<> X, doubles y, doubles_matrix<> U, doubles w, int nthreads)
{
    // We compute the following
    // - (UX)'(UX)
    // - (UX)'y

    int N = U.nrow();
    int K1 = U.ncol();

    bool isX = X.nrow() > 1;
    int K2 = isX ? X.ncol() : 0;

    bool isWeight = w.size() > 1;

    writable::doubles_matrix<> UXtUX(K2 + K1, K2 + K1);
    writable::doubles UXty(K2 + K1);

    writable::doubles_matrix<> wU(U);
    if (isWeight)
    {
        for (int k = 0; k < K1; ++k)
        {
            for (int i = 0; i < N; ++i)
            {
                wU(i, k) *= w[i];
            }
        }
    }

    writable::list res;

    // (UX)'y
    for (int k = 0; k < K2; ++k)
    {
        UXty[k + K1] = Xty[k];
    }

    // TODO: OMP functions
    // #pragma omp parallel for num_threads(nthreads)
    for (int k = 0; k < K1; ++k)
    {
        double val = 0;
        for (int i = 0; i < N; ++i)
        {
            val += y[i] * wU(i, k);
        }

        UXty[k] = val;
    }

    res["UXty"] = UXty;

    // (UX)'(UX)

    if (sparse_check(X) == false)
    {
        mp_ZXtZX(UXtUX, XtX, X, U, wU, nthreads);
        res["UXtUX"] = UXtUX;
    }
    else
    {
        vector<int> n_j(K2 + !isX, 0);
        vector<int> start_j(K2 + !isX + 1, 0);
        vector<int> all_i;
        vector<double> x;

        set_sparse(n_j, start_j, all_i, x, X, w);

        mp_sparse_ZXtZX(UXtUX, XtX, n_j, start_j, all_i, x, X, U, wU, nthreads);
        res["UXtUX"] = UXtUX;
    }

    return res;
}

[[cpp11::register]] doubles cpp_iv_resid(doubles resid_2nd, doubles coef, SEXP resid_1st, bool is_int, int nthreads)
{

    int N = resid_2nd.size();
    int K = Rf_length(resid_1st);

    writable::doubles iv_resid(resid_2nd);
    vector<int> bounds = set_parallel_scheme(N, nthreads);

    if (K == 1)
    {
        double *p_r = REAL(VECTOR_ELT(resid_1st, 0));

        // TODO: OMP functions
        // #pragma omp parallel for num_threads(nthreads)
        for (int t = 0; t < nthreads; ++t)
        {
            for (int i = bounds[t]; i < bounds[t + 1]; ++i)
            {
                iv_resid[i] -= coef[0 + is_int] * p_r[i];
            }
        }
    }
    else
    {

        vector<double *> p_p_r(K);
        for (int k = 0; k < K; ++k)
        {
            p_p_r[k] = REAL(VECTOR_ELT(resid_1st, k));
        }

        // TODO: OMP functions
        // #pragma omp parallel for num_threads(nthreads)
        for (int t = 0; t < nthreads; ++t)
        {
            for (int k = 0; k < K; ++k)
            {
                double *p_r = p_p_r[k];
                for (int i = bounds[t]; i < bounds[t + 1]; ++i)
                {
                    iv_resid[i] -= coef[k + is_int] * p_r[i];
                }
            }
        }
    }

    return iv_resid;
}
