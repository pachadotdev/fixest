/*
 _ __ ___   _  ___   ___
| '_ ` _ \ | |/ __| / __|
| | | | | || |\__ \| (__
|_| |_| |_||_||___/ \___|

Original Author: Laurent R. Berge
Refactored by Mauricio "Pacha" Vargas Sepulveda starting in Jun 2022

Three groups of non-parallel functions:
1) simple functions that are faster than base R (because more specific)
2) function to obtain the FE coefficients after the estimation is done
3) functions to lag variables
*/

// TODO: OMP functions
#pragma once

#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <vector>
#include <functional>

#include <stdint.h>

#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

using namespace cpp11;
using std::vector;

class simple_vec_double
{
    simple_vec_double() = delete;
    double *px_double = nullptr;
    int *px_int = nullptr;
    int n;
    bool is_real;

public:
    simple_vec_double(SEXP x);
    double operator[](int);
};
