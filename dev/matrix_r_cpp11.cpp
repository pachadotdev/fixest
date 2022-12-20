#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/list.hpp>
#include <cpp11/matrix.hpp>
#include <cpp11/logicals.hpp>

using namespace cpp11;

[[cpp11::register]] doubles_matrix<> matrix_r_cpp11(doubles_matrix<> X) {
    int K = X.ncol();
    writable::doubles_matrix<> R(K, K);
    double min_norm = X(0, 0);

    for (int j = 0; j < K; ++j) {
        double R_jj = X(j, j);
        for (int k = 0; k < j; ++k) {
            R_jj -= R(k, j) * R(k, j);
        }

        if (min_norm > R_jj) { min_norm = R_jj; }
        R_jj = sqrt(R_jj);
        R(j, j) = R_jj;

        for (int i = j + 1; i < K; ++i) {
            double value = X(i, j);
            for (int k = 0; k < j; ++k) {
                value -= R(k, i) * R(k, j);
            }
            R(j, i) = value / R_jj;
        }
    }

    return R;
}
