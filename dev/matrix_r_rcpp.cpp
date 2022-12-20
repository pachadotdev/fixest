#include <Rcpp.h>
#include <cmath>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix matrix_r_rcpp(NumericMatrix X){
    int K = X.ncol();
    NumericMatrix R(K, K);
    double min_norm = X(0, 0);

    for(int j=0 ; j<K ; ++j){
        double R_jj = X(j, j);
        for(int k=0 ; k<j ; ++k) {
            R_jj -= R(k, j) * R(k, j);
        }

        if(min_norm > R_jj) min_norm = R_jj;
        R_jj = sqrt(R_jj);
        R(j, j) = R_jj;

        for(int i=j+1 ; i<K ; ++i) {
            double value = X(i, j);
            for(int k=0 ; k<j ; ++k) {
                value -= R(k, i) * R(k, j);
            }
            R(j, i) = value/R_jj;
        }
    }

    return R;
}
