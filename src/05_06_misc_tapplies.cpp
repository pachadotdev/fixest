#include "05_0_misc.h"

// Function to get the conditional sum of a matrix
[[cpp11::register]] doubles_matrix<> cpp_tapply_sum_(int Q, doubles_matrix<> x, integers dum)
{
    // Q: nber of classes
    // N: nber of observations
    // x: a matrix
    // dum: the N vector of clusters

    int N = x.nrow();
    int K = x.ncol();

    writable::doubles_matrix<> res(Q, K);

    // fill with 0
    for (int q = 0; q < Q; q++)
    {
        for (int k = 0; k < K; k++)
        {
            res(q, k) = 0;
        }
    }

    int i, q, k;

    for (i = 0; i < N; i++)
    {
        q = dum[i] - 1; // we take 1 off => different indexation in C

        for (k = 0; k < K; k++)
        {
            res(q, k) += x(i, k);
        }
    }

    return (res);
}

// Function to get the conditional sum of a vector
[[cpp11::register]] doubles cpp_tapply_vsum_(int Q, doubles x, integers dum)
{
    // Q: nber of classes
    // x: a matrix
    // dum: the N vector of clusters

    int N = x.size();

    writable::doubles res(Q);
    int i, q;

    for (i = 0; i < N; i++)
    {
        q = dum[i] - 1; // we take 1 off => different indexation in C
        res[q] += x[i];
    }

    return (res);
}
