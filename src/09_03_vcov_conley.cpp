#include "09_0_vcov.h"

mat_row_scheme::mat_row_scheme(doubles_matrix<> &x)
{

    this->N = x.nrow();
    this->K = x.ncol();
    n_total = N * K;

    // N and K are int64, so product is OK
    mat.resize(n_total);

    // filling scheme is slow, memcpy cannot be leveraged
    for (int64_t i = 0; i < N; ++i)
    {
        for (int64_t k = 0; k < K; ++k)
        {
            mat[i * K + k] = x(i, k);
        }
    }
}

mat_row_scheme::mat_row_scheme(mat_row_scheme &x)
{
    // copy

    this->N = x.nrow();
    this->K = x.ncol();

    // N and K are int64, so product is OK
    n_total = N * K;
    mat.resize(n_total);

    // std::copy can be used => fast
    std::copy(x.mat.begin(), x.mat.end(), mat.begin());
}

mat_row_scheme::mat_row_scheme(int N_in, int K_in)
{

    this->N = N_in;
    this->K = K_in;
    n_total = N * K;

    // N and K are int64, so product is OK
    mat.resize(n_total);

    std::fill(mat.begin(), mat.end(), 0);
}

// S: scores
// lon_rad/lat_rad: longitude/latitude
//  => the lat_rad and lon_rad **must** be in radians!!!!
//
// IMPORTANT: S, lon_rad, lat_rad must be sorted by lat
//
// cutoff, in km
// kernel type: uniform only, it does not matter
// distance:
// - 1: spherical
// - 2: triangular

// Notes to self:
// [border problem]
// - when we use latitude as main dimension, there are no border problems
// - border problems (that -170 is close to 170) happen only for longitude
// - using only the latitude is the way to go, even though the data is usually more scattered across longitude
// and it would save time to use longitude as the main cutoff, that would imply just too many computations
//
// [longitude cutoff]
// - the longitude cutoff depends on the latitude of the points
// - I initially used only the latitude of i
// - but there are two points with different latitudes => thus that breaks symmetry!
// - in the end I took the avg latitudes
// - I think that in some rare instances there can still be cases where:
//   + dist_ij < cutoff, but..
//   + dist_lon_rad_ij > cutoff_lon_mean_ij (and possibly: dist_lon_rad_ij > cutoff_lon_i OR dist_lon_rad_ij > cutoff_lon_j)
// - but it should be super edge cases, at the very limit... so not really important
//
// [using longitude]
// - given the problems mentioned above, I should NEVER use longitude as the main dimension
// - Laurent, remember you first started with longitude in main and although you can't mention here
// all the problems that you encountered, you did encounter them. So don't.
[[cpp11::register]] doubles_matrix<> cpp_vcov_conley(doubles_matrix<> S, doubles lon_rad, doubles lat_rad,
                                                     const int distance, const double cutoff, int nthreads)
{
    if (distance <= 0 || distance > 2)
    {
        stop("'distance' is not valid (internal error).");
    }

    int K = S.ncol();
    int N = S.nrow();

    mat_row_scheme scores(S);

    // utilities
    writable::doubles cos_lat(N);
    for (int i = 0; i < N; ++i)
    {
        cos_lat[i] = cos(lat_rad[i]);
    }

    mat_row_scheme cum_scores(scores);

    cum_scores.scale(0.5);

    // 1 lat degree is approx 111km
    const double lat_cutoff_rad = degree_to_radian(cutoff / 111);
    const double lon_cutoff_rad_factor = degree_to_radian(cutoff / 111);

    // cutoff_rad_sq used when distance == 2
    const double cutoff_rad_sq = to_sq(degree_to_radian(cutoff) / 111);

    // TODO: OMP functions
    // #pragma omp parallel for num_threads(nthreads)
    for (int i = 0; i < N; ++i)
    {

        const double lon_rad_i = lon_rad[i];
        const double lat_rad_i = lat_rad[i];
        const double cos_lat_i = cos_lat[i];

        // 1 lon degree is approx cos(lat)*111km
        double lon_cutoff_rad = 0;

        bool ok = false;
        double dist_lon_rad = 0;
        double dist_lat_rad = 0;
        double dist = 0;
        double cos_lat_mean = 0;

        for (int i2 = i + 1; i2 < N; ++i2)
        {
            dist_lat_rad = fabs_lat(lat_rad[i2], lat_rad_i);
            if (dist_lat_rad > lat_cutoff_rad)
            {
                break;
            }

            dist_lon_rad = fabs_lon(lon_rad[i2], lon_rad_i);

            // we take the avg lat for symmetry
            cos_lat_mean = cos((lat_rad_i + lat_rad[i2]) / 2);
            lon_cutoff_rad = lon_cutoff_rad_factor / cos_lat_mean;

            if (dist_lon_rad > lon_cutoff_rad)
            {
                continue;
            }

            if (distance == 1)
            {
                dist = dist_km(lon_rad_i, lat_rad_i, cos_lat_i,
                               lon_rad[i2], lat_rad[i2], cos_lat[i2]);
                ok = dist <= cutoff;
            }
            else if (distance == 2)
            {
                // dist: in radian, and put to square
                dist = to_sq(dist_lat_rad) + to_sq(cos_lat_mean * dist_lon_rad);
                ok = dist <= cutoff_rad_sq;
            }

            // if(cpt++ < 10) Rcout << "i = " << i << ", j = " << i2 << ", dij_lon = " << dist_lon_rad << ", dij_lat = " << dist_lat_rad << ", d_ij = " << dist << ", w = " << weight << "\n";

            if (ok)
            {
                // ++n_done;
                for (int k = 0; k < K; ++k)
                {
                    cum_scores(i, k) += scores(i2, k);
                }
            }
        }
    }

    // Rcout << "scores:\n";
    // scores.check();
    // Rcout << "cum scores:\n";
    // cum_scores.check();

    // Now the matrix multiplications

    // Rcout << "total done: " << ++n_done << "\n";

    writable::doubles_matrix<> res(K, K);

    for (int i = 0; i < N; ++i)
    {
        for (int k1 = 0; k1 < K; ++k1)
        {
            for (int k2 = 0; k2 < K; ++k2)
            {
                res(k1, k2) += scores(i, k1) * cum_scores(i, k2);
            }
        }
    }

    // we add the transpose
    double tmp = 0;
    for (int k1 = 0; k1 < K; ++k1)
    {
        for (int k2 = k1; k2 < K; ++k2)
        {
            if (k1 == k2)
            {
                res(k1, k2) *= 2;
            }
            else
            {
                tmp = res(k1, k2);
                res(k1, k2) += res(k2, k1);
                res(k2, k1) += tmp;
            }
        }
    }

    return res;
}
