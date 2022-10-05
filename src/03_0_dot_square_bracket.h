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

// TODO: OMP functions #pragma once

#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <vector>
#include <cstring>

#include <stdint.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

using namespace cpp11;
using std::vector;
using std::string;
using std::strlen;

vector<int> set_parallel_scheme(int N, int nthreads);

inline bool is_dsb_open(const char *str, int i, int n);
