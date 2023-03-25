/*
                __
               / _|
  __ _  _   _ | |_
 / _` || | | ||  _|
| (_| || |_| || |
 \__, | \__,_||_|
    | |
    |_|

Original Author: Laurent R. Berge
Refactored by Mauricio "Pacha" Vargas Sepulveda starting in Jun 2022

QUF stands for "quick unclass factor", an operation consisting
in transforming a vector of arbitrary values into an integer
vector ranging from 1 to the number of unique values (i.e:
  a)  (50, 55, 32, 12, 12) => (3, 4, 2, 1, 1)
  b)  ("a", "d", "a", "b") => (1, 3, 1, 2)

 The code here is used to quf vectors of integers, floats
 or strings. I convert any other type of identifier to
 character if they're not numeric before getting into this
 function.

 All we care in qufing is to get a vector from 1 to n the number
 of unique values. We don't care about the order (e.g., in
 example a) the result (1,2,3,4,4) would have been fine).

 Therefore, qufing integer vectors of small range is both
 very simple and very fast: we just have to set a lookup table.
 For integers of large range, setting up a large
 lookup table can be too costly, therefore I sort the vector
 first and then create the new values.

 The process is similar for doubles: first sort, then get the
 new values.

 Finally for string vectors, I use R's feature of stocking
 character strings in a unique location. They are therefore
 uniquely identified by their pointer. I then transform the
 pointer to int, and sort accordingly.

 The sorting, when sorting is required, is radix-based. My code
 has greatly benefited from two sources which clarified a lot
 of implicit things:
 - http://codercorner.com/RadixSortRevisited.htm (Pierre Terdiman)
 - http://stereopsis.com/radix.html (Michael Herf)

 As said before, all that matters for qufing is that identical
 values are consecutive after the sort. Thus, the endianness of
 the system is of no importance: whatever the order of bytes
 on which we sort, we obtain what we want.
*/

#pragma once

#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <vector>
#include <numeric>

#include <stdint.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

using namespace cpp11;
using std::vector;
using std::string;
using std::uintptr_t;
using std::fill;
using std::accumulate;

void quf_refactor(int *px_in, int x_size, integers &obs2keep, int n, int *x_uf, vector<double> &x_unik, vector<int> &x_table);

void quf_single(void *px_in, string &x_type, int n, int *x_uf, vector<double> &x_unik);

void quf_refactor_table_sum_single(int n, int *quf_old, int *quf_new, vector<bool> &obs_removed,
                                   vector<double> &x_unik, vector<double> &x_unik_new, vector<double> &x_removed,
                                   vector<int> &x_table, double *py, vector<double> &sum_y, bool do_sum_y,
                                   bool rm_1, bool rm_single, vector<bool> &id_pblm, bool check_pblm,
                                   bool *pstop_now, int n_FE);
