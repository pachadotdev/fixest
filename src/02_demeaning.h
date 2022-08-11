/*
     _                                       _               
    | |                                     (_)              
  __| |  ___  _ __ ___    ___   __ _  _ __   _  _ __    __ _ 
 / _` | / _ \| '_ ` _ \  / _ \ / _` || '_ \ | || '_ \  / _` |
| (_| ||  __/| | | | | ||  __/| (_| || | | || || | | || (_| |
 \__,_| \___||_| |_| |_| \___| \__,_||_| |_||_||_| |_| \__, |
                                                        __/ |
                                                       |___/ 

Original Author: Laurent R. Berge
Refactored by Mauricio "Pacha" Vargas Sepulveda starting in Jun 2022

Workhorse for feols and feglm.

It demeans any variable given in input, the algortihm is not
a real demeaning algorithm.
It is in fact identical as obtaining the optimal set of
cluster coefficients (i.e. fixed-effects) in a ML context with
a Gaussian likelihood (as described in Berge, 2018).

This way we can leverage the powerful Irons and Tuck acceleration
algorithm. For simple cases, this doesn't matter so much.
But for messy data (employee-company for instance), it matters
quite a bit.

In terms of functionality it accommodate weights and coefficients
with varying slopes.

Of course any input is **strongly** checked before getting into
this function.

I had to apply a trick to accommodate user interrupt in a
parallel setup. It costs a bit, but it's clearly worth it.
*/

/* 
CONVENTIONS

Suffixes
_Q: the vector is of length Q
_I (for Identifier): the vector is of length the total number of FE identifiers
_T: a scalar representing the Total of something, usually the sum of a _Q object
_C: the vector is of length the number of coefficients
VS_C: the vector is of length the number of coefficients for the varying slopes
noVS_C: the vector is of length the number of coefficients for regular fixed-effects (no VS!)
_N: the vector is of length n_obs

PREFIXES
p_ means a pointer
nb_id: refers to the fixed-effects IDs
nb_coef: refers to the fixed-effects coefficients. Coef and id are identical in **absence** of varying slopes.
vs: means varying slopes
 */

// TODO: Next version => clean c++ code, use only sMat

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
using std::fabs;

// COMMON FUNS

bool continue_criterion(double a, double b, double diffMax);

bool stopping_criterion(double a, double b, double diffMax);

bool update_X_IronsTuck(int nb_coef_no_K, vector<double> &X,
						const vector<double> &GX, const vector<double> &GGX,
						vector<double> &delta_GX, vector<double> &delta2_X);

// DEMEANING

class sVec
{
    double *p_dble = nullptr;
    int *p_int = nullptr;

public:
    // several constructors

    // is_int public member
    bool is_int = false;

    sVec(){};
    sVec(SEXP);
    sVec(double *p_x) : p_dble(p_x), is_int(false){};
    sVec(int *p_x) : p_int(p_x), is_int(true){};
    sVec(std::nullptr_t){};

    double operator[](int i)
    {
        if (is_int)
            return static_cast<double>(p_int[i]);
        return p_dble[i];
    }
};

class sMat
{

    std::vector<sVec> p_sVec;
    int n = 0;
    int K = 0;

    sMat() = delete;

public:
    sMat(SEXP);

    int nrow() { return n; };
    int ncol() { return K; };

    sVec operator[](int);
    double operator()(int, int);
};

class FEClass
{

    int Q;
    int n_obs;
    int slope_index;
    bool is_weight;
    bool is_slope;

    // Dense vectors that we populate in the class, and their associated pointers
    vector<double> eq_systems_VS_C;
    vector<double *> p_eq_systems_VS_C;

    vector<double> sum_weights_noVS_C;
    vector<double *> p_sum_weights_noVS_C;

    // p_fe_id: pointers to the fe_id vectors
    // p_vs_vars: pointers to the VS variables
    // p_weights: pointer to the weight vector
    // eq_systems_VS_C: vector stacking all the systems of equations (each system is of size n_coef * n_vs * n_vs)
    // p_eq_systems_VS_C: pointer to the right equation system. Of length Q.
    vector<int *> p_fe_id;
    vector<sVec> p_vs_vars;
    double *p_weights = nullptr;

    vector<bool> is_slope_Q;
    vector<bool> is_slope_fe_Q;

    vector<int> nb_vs_Q;
    vector<int> nb_vs_noFE_Q;

    int *nb_id_Q;

    vector<int> coef_start_Q;

    // internal functions
    void compute_fe_coef_internal(int, double *, bool, sVec, double *, double *);
    void compute_fe_coef_2_internal(double *, double *, double *, bool);
    void add_wfe_coef_to_mu_internal(int, double *, double *, bool);

    public:
    // Utility class: Facilitates the access to the VS variables
    class simple_mat_of_vs_vars
    {
        int K_fe;
        vector<sVec> pvars;

    public:
        simple_mat_of_vs_vars(const FEClass *, int);
        double operator()(int, int);
    };

    int nb_coef_T;
    vector<int> nb_coef_Q;

    // constructor:
    FEClass(int n_obs, int Q, SEXP r_weights, SEXP fe_id_list, SEXP r_nb_id_Q, SEXP table_id_I, SEXP slope_flag_Q, SEXP slope_vars_list);

    // functions
    void compute_fe_coef(double *fe_coef, sVec &mu_in_N);
    void compute_fe_coef(int q, double *fe_coef, double *sum_other_coef_N, double *in_out_C);

    void add_wfe_coef_to_mu(int q, double *fe_coef_C, double *out_N);
    void add_fe_coef_to_mu(int q, double *fe_coef_C, double *out_N);

    void compute_fe_coef_2(double *fe_coef_in_C, double *fe_coef_out_C, double *fe_coef_tmp, double *in_out_C);

    void add_2_fe_coef_to_mu(double *fe_coef_a, double *fe_coef_b, double *in_out_C, double *out_N, bool update_beta);

    void compute_in_out(int q, double *in_out_C, sVec &in_N, double *out_N);
};

class simple_mat_with_id
{
    // => Access to one of the n_coef matrices of size n_vs x n_vs; all stacked in a single vector
    //
    // px0: origin of the vector (which is of length n_coef * n_vs * n_vs)
    // px_current: working location of the n_vs x n_vs matrix
    // n: n_vs
    // n2: explicit
    // id_current: current identifier. The identifiers range from 0 to (n_coef - 1)

    simple_mat_with_id() = delete;

    double *px0;
    double *px_current;
    int nrow, ncol, n_total, id_current = 0;

public:
    simple_mat_with_id(double *px_in, int nrow_in) : px0(px_in), px_current(px_in), nrow(nrow_in), ncol(nrow_in), n_total(nrow * ncol){};
    simple_mat_with_id(double *px_in, int nrow_in, int ncol_in) : px0(px_in), px_current(px_in), nrow(nrow_in), ncol(ncol_in), n_total(nrow * ncol){};
    double &operator()(int id, int i, int j);
    double &operator()(int id, int i);
};

struct PARAM_DEMEAN
{
    int n_obs;
    int Q;
    int nb_coef_T;
    int iterMax;
    double diffMax;

    // iterations
    int *p_iterations_all;

    // vectors of pointers
    vector<sVec> p_input;
    vector<double *> p_output;

    // saving the fixed effects
    bool save_fixef;
    double *fixef_values;

    // FE information
    FEClass *p_FE_info;

    // stopflag
    bool *stopnow;
    int *jobdone;
};

void compute_fe_gnl(double *p_fe_coef_origin, double *p_fe_coef_destination,
                    double *p_sum_other_means, double *p_sum_in_out, PARAM_DEMEAN *args);

void stayIdleCheckingInterrupt(bool *stopnow, vector<int> &jobdone, int n_vars, int *counterInside);
