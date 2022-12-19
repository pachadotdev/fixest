/**************************************************************
 * ___________________                                         *
 * || OLS Functions ||                                         *
 * -------------------                                         *
 *                                                             *
 * Author: Laurent R. Berge                                    *
 *                                                             *
 * Set of functions to perform OLS estimations.                *
 *                                                             *
 * In general, the functions are slower than BLAS ones         *
 * but they have the merit of doing exacly what I want.        *
 *                                                             *
 **************************************************************/

#include <Rcpp.h>
#include <cmath>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

void invert_tri(NumericMatrix &R, int K, int nthreads = 1){

    // Startegy: we invert by bands (b) => better for parallelization

    // initialization of R prime
    for(int i=0 ; i<K ; ++i){
        for(int j=i+1 ; j<K ; ++j){
            R(j, i) = R(i, j);
        }
    }

    // b0
    for(int i=0 ; i<K ; ++i){
        R(i, i) = 1/R(i, i);
    }

    // Check for interrupts
    // number of computations is (K - b) * (b + 1) => max is (K + 1)**2 / 2
    double flop = (K + 1) * (K + 1) / 2.0;
    int iterSecond = ceil(2000000000 / flop / 5); // nber iter per 1/5 second

    for(int b=1 ; b<K ; ++b){

        if(b % iterSecond == 0){
            // Rprintf("check\n");
            R_CheckUserInterrupt();
        }

#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
        for(int i=0 ; i<K-b ; ++i){
            double numerator = 0;
            int col = i+b;
            for(int k=i+1 ; k<=col ; ++k){
                numerator -= R(k, i) * R(k, col);
            }
            R(i, col) = numerator * R(i, i);
        }
    }

}

void tproduct_tri(NumericMatrix &RRt, NumericMatrix &R, int nthreads = 1){

    int K = RRt.ncol();

    // initialization of R prime
    for(int i=0 ; i<K ; ++i){
        for(int j=i+1 ; j<K ; ++j){
            R(j, i) = R(i, j);
        }
    }

    // Check for interrupts
    // we do the same as for the invert_tri
    double flop = (K + 1) * (K + 1) / 2.0;
    int iterSecond = ceil(2000000000 / flop / 5); // nber iter per 1/5 second
    int n_iter_main = 0;

#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
    for(int i=0 ; i<K ; ++i){

        if(omp_get_thread_num() == 0 && n_iter_main % iterSecond == 0){
            // Rprintf("check\n");
            R_CheckUserInterrupt();
            ++n_iter_main;
        }

        for(int j=i ; j<K ; ++j){

            double value = 0;
            int k_start = i < j ? j : i;
            for(int k=k_start ; k<K ; ++k){
                value += R(k, j) * R(k, i);
            }
            RRt(i, j) = value;
            RRt(j, i) = value;
        }
    }

}


// [[Rcpp::export]]
List cpp_cholesky(NumericMatrix X, double tol = 1.0/100000.0/100000.0, int nthreads = 1){
    // X est symetrique, semi definie positive
    // rank-revealing on-the-fly

    List res;

    int K = X.ncol();

    NumericMatrix R(K, K);
    LogicalVector id_excl(K);
    int n_excl = 0;

    // we check for interrupt every 1s when it's the most computationally intensive
    // at each iteration we have K * (j+1) - j**2 - 2*j - 1 multiplications
    // max => K**2/4
    double flop = K * K / 4.0;
    int iterSecond = ceil(2000000000 / flop / 5); // nber iter per 1/5 second
    double min_norm = X(0, 0);

    for(int j=0 ; j<K ; ++j){

        if(j % iterSecond == 0){
            // Rprintf("check\n");
            R_CheckUserInterrupt();
        }

        // implicit pivoting (it's like 0-rank variables are stacked in the end)

        double R_jj = X(j, j);
        for(int k=0 ; k<j ; ++k){

            if(id_excl[k]) continue;

            R_jj -= R(k, j) * R(k, j);
        }

        if(R_jj < tol){
            n_excl++;
            id_excl[j] = true;

            // Rcout << "excluded: " << j << ", L_jj: " << L_jj << "\n";

            // Corner case, may happen:
            if(n_excl == K){
                List res;
                res["all_removed"] = true;
                return res;
            }

            continue;
        }

        if(min_norm > R_jj) min_norm = R_jj;

        R_jj = sqrt(R_jj);

        R(j, j) = R_jj;

#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
        for(int i=j+1 ; i<K ; ++i){

            double value = X(i, j);
            for(int k=0 ; k<j ; ++k){
                if(id_excl[k]) continue;

                value -= R(k, i) * R(k, j);
            }
            R(j, i) = value/R_jj;
        }
    }

    res["R"] = R;

    // we reconstruct the R matrix, otherwise it's just a mess to handle excluded variables on the fly
    // changes are in place

    if(n_excl > 0){
        int n_j_excl = 0;

        // The first chunk of the matrix is identical
        // we find when it starts being different
        int j_start = 0;
        for( ; !id_excl[j_start] ; ++j_start);

        for(int j=j_start ; j<K ; ++j){

            if(id_excl[j]){
                ++n_j_excl;
                continue;
            }

            int n_i_excl = 0;
            for(int i=0 ; i<=j ; ++i){
                if(id_excl[i]){
                    ++n_i_excl;
                    continue;
                }
                R(i - n_i_excl, j - n_j_excl) = R(i, j);
            }
        }

        K -= n_excl;
    }

    res["R"] = R;
    res["K"] = K;

    // Inversion of R, in place
    invert_tri(R, K, nthreads);

    NumericMatrix XtX_inv(K, K);
    tproduct_tri(XtX_inv, R, nthreads);

    res["iR"] = R;
    res["XtX_inv"] = XtX_inv;
    res["id_excl"] = id_excl;
    res["min_norm"] = min_norm;

    return res;
}
