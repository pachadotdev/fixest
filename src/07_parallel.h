/*******************************************************************
 * ________________________                                        *
 * || Parallel functions ||                                        *
 * ------------------------                                        *
 *                                                                 *
 * Author: Laurent R. Berge                                        *
 *                                                                 *
 * Group of functions doing simple things... but in parallel.      *
 *                                                                 *
 * The functions don't do much more than what their names suggest. *
 *                                                                 *
 ******************************************************************/

#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <cfloat>
#include <math.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#include <pthread.h>
#else
#define omp_get_max_threads() 0
#endif
#include <cmath>
#include <stdio.h>
#include <Rmath.h>

using namespace cpp11;

// This file contains misc fixest functions parallelized with the omp library

// Regarding fork detection => I don't know enough yet
// The code below doesn't seem to work. And I can't check since I'm on Windows....
// Safer not to include it for now.

// static bool fixest_in_fork = false;
//
// // Trick taken from data.table to detect forking
// // Actually I don't know if it works, and I can't check...
// void when_fork() {
//   fixest_in_fork = true;
// }
//
// void after_fork() {
//   fixest_in_fork = false;
// }
