/*
     _       _
    | |     | |
  __| | ___ | |__
 / _` |/ __|| '_ \
| (_| |\__ \| |_) |
 \__,_||___/|_.__/

Original Author: Laurent R. Berge
Refactored by Mauricio "Pacha" Vargas Sepulveda starting in Jun 2022
*/

#include "00_common.hpp"

#pragma once

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

vector<int> set_parallel_scheme(int N, int nthreads);

inline bool is_dsb_open(const char *str, int i, int n);
