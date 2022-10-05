#include "05_0_misc.h"

[[cpp11::register]] doubles cpp_lgamma(doubles x)
{
    // simple function to compute lgamma of a vector

    int n = x.size();
    writable::doubles res(n);

    for (int i = 0; i < n; i++)
    {
        res[i] = lgamma(x[i]);
    }

    return (res);
}

[[cpp11::register]] doubles cpp_log_a_exp(double a, doubles mu, doubles exp_mu)
{
    // faster this way

    int n = mu.size();
    writable::doubles res(n);

    for (int i = 0; i < n; i++)
    {
        if (mu[i] < 200)
        {
            res[i] = log(a + exp_mu[i]);
        }
        else
        {
            res[i] = mu[i];
        }
    }

    return (res);
}
