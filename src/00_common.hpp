// #include <cstring>
// #include <functional>
// #include <algorithm>
#include <Rmath.h>
#include <float.h>
#include <math.h>
#include <stdint.h>

#include <cpp11.hpp>
#include <iostream>
#include <numeric>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

using namespace cpp11;

bool continue_criterion(double a, double b, double diffMax);

bool stopping_criterion(double a, double b, double diffMax);

bool update_X_IronsTuck(int nb_coef_no_K, std::vector<double> &X,
                        const std::vector<double> &GX,
                        const std::vector<double> &GGX,
                        std::vector<double> &delta_GX,
                        std::vector<double> &delta2_X);
