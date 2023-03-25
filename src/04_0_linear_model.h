/*
 _
| |
| | _ __ ___
| || '_ ` _ \
| || | | | | |
|_||_| |_| |_|

Original Author: Laurent R. Berge
Refactored by Mauricio "Pacha" Vargas Sepulveda starting in Jun 2022

Set of functions to perform OLS estimations.

In general, the functions are slower than BLAS ones
but they have the merit of doing exacly what I want.
*/

#pragma once

#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <vector>

#include <stdint.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

using namespace cpp11;
using std::vector;

vector<int> set_parallel_scheme(int N, int nthreads);

bool sparse_check(const doubles_matrix<> &X);

void set_sparse(std::vector<int> &n_j, std::vector<int> &start_j, std::vector<int> &all_i, std::vector<double> &x, const doubles_matrix<> &X, const doubles &w);

void mp_XtX(writable::doubles_matrix<> &XtX, const doubles_matrix<> &X, const doubles_matrix<> &wX, int nthreads);

void mp_Xty(writable::doubles &Xty, const doubles_matrix<> &X, const double *y, int nthreads);

void mp_ZXtu(writable::doubles &ZXtu, const doubles_matrix<> &X, const doubles_matrix<> &Z, const double *u, int nthreads);

void mp_ZXtZX(writable::doubles_matrix<> &ZXtZX, const doubles_matrix<> &XtX, const doubles_matrix<> &X, const doubles_matrix<> &Z, const doubles_matrix<> &wZ, int nthreads);

void mp_sparse_XtX(writable::doubles_matrix<> &XtX, const std::vector<int> &n_j, const std::vector<int> &start_j, const std::vector<int> &all_i, const std::vector<double> &x, const doubles_matrix<> &X, int nthreads);

void mp_sparse_Xty(writable::doubles &Xty, const std::vector<int> &start_j, const std::vector<int> &all_i, const std::vector<double> &x, const double *y, int nthreads);

void mp_sparse_ZXtu(writable::doubles &ZXtu, const std::vector<int> &start_j, const std::vector<int> &all_i, const std::vector<double> &x, const double *u, const doubles_matrix<> &X, const doubles_matrix<> &wZ, int nthreads);

void mp_sparse_ZXtZX(writable::doubles_matrix<> &ZXtZX, const doubles_matrix<> &XtX, const std::vector<int> &n_j, const std::vector<int> &start_j, const std::vector<int> &all_i, const std::vector<double> &x, const doubles_matrix<> &X, const doubles_matrix<> &Z, const doubles_matrix<> &wZ, int nthreads);
