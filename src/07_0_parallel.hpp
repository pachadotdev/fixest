/*
                            _  _        _
                           | || |      | |
 _ __    __ _  _ __   __ _ | || |  ___ | |
| '_ \  / _` || '__| / _` || || | / _ \| |
| |_) || (_| || |   | (_| || || ||  __/| |
| .__/  \__,_||_|    \__,_||_||_| \___||_|
| |
|_|

Original Author: Laurent R. Berge
Refactored by Mauricio "Pacha" Vargas Sepulveda starting in Jun 2022

This file contains misc fixest functions parallelized with the omp library

Regarding fork detection => I don't know enough yet
The code below doesn't seem to work. And I can't check since I'm on Windows...
Safer not to include it for now.

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
*/

#include "00_common.hpp"

#pragma once

std::vector<int> set_parallel_scheme_bis(int N, int nthreads);
