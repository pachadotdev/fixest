/***********************************************************************************
 * ________________                                                                *
 * || Misc. funs ||                                                                *
 * ----------------                                                                *
 *                                                                                 *
 * Author: Laurent R. Berge                                                        *
 *                                                                                 *
 * Three groups of non-parallel functions:                                         *
 *   1) simple functions that are faster than base R (because more specific)       *
 *   2) function to obtain the FE coefficients after the estimation is done        *
 *   3) functions to lag variables                                                 *
 *                                                                                 *
 * Nothing much to say, everything is pretty explicit.                             *
 *                                                                                 *
 **********************************************************************************/

#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <functional>

using namespace cpp11;
using namespace std;

[[cpp11::register]]
doubles cpp_lgamma(doubles x){
    // simple function to compute lgamma of a vector

    int n = x.size();
    writable::doubles res(n);

    for(int i=0 ; i<n ; i++){
        res[i] = lgamma(x[i]);
    }

    return(res);
}

[[cpp11::register]]
doubles cpp_log_a_exp(double a, doubles mu, doubles exp_mu){
    // faster this way

    int n = mu.size();
    writable::doubles res(n);

    for(int i=0 ; i<n ; i++){
        if(mu[i] < 200){
            res[i] = log(a + exp_mu[i]);
        } else {
            res[i] = mu[i];
        }
    }

    return(res);
}

[[cpp11::register]]
doubles cpp_partialDerivative_other(int iterMax,
                                    int Q,
                                    int N,
                                    double epsDeriv,
                                    doubles ll_d2,
                                    doubles dx_dother,
                                    doubles init,
                                    integers_matrix<> dumMat,
                                    integers nbCluster){
    // takes in:
    // dumMat: the matrix of dummies (n X c) each obs => cluster // must be in cpp index!!!
    // init: the initialisation of the sum of derivatives vector
    // ll_d2: the second derivative
    // dx_dother: the vector of dx_dother

    int iter;

    int i, q, c;
    int index;
    int sum_cases=0;
    bool ok;
    double new_value;

    // OLD: int start[Q], end[Q];
    // NEW: using ISO C++
    std::vector <int> start, end;
    start.reserve(Q);
    end.reserve(Q);

    for(q=0 ; q<Q ; q++) {
        // the total number of clusters (eg if man/woman and 10 countries: total of 12 cases)
        sum_cases += nbCluster[q];
        if(q == 0){
            start[q] = 0;
            end[q] = nbCluster[q];
        } else {
            start[q] = start[q-1] + nbCluster[q-1];
            end[q] = end[q-1] + nbCluster[q];
        }
    }

    writable::doubles clusterDeriv(sum_cases); // the derivatives per cluster, stacked in one vector
    writable::doubles sum_lld2(sum_cases);

    // Creation of the sum_lld2
    for (i=0 ; i<N ; i++) {
        for (q=0 ; q<Q ; q++) {
            index = start[q] + dumMat(i, q);
            sum_lld2[index] += ll_d2[i];
        }
    }

    // the result to return
    writable::doubles S;
    for(i=0 ; i<N ; i++){
        S[i] = init[i];
    }

    ok = true;
    iter = 0;
    while( ok & (iter<iterMax) ){
        iter++;
        ok = false;

        for(q=0 ; q<Q ; q++){
            R_CheckUserInterrupt();

            // init of the vector
            for(c=start[q] ; c<end[q] ; c++){
                clusterDeriv[c] = 0;
            }

            for(i=0 ; i<N ; i++){
                index = start[q] + dumMat(i, q);
                clusterDeriv[index] += dx_dother[i] + S[i]*ll_d2[i];
            }

            // on finit de calculer clusterDeriv + controle
            for(c=start[q] ; c<end[q] ; c++){
                new_value = (-1.0) * clusterDeriv[c] / sum_lld2[c];
                clusterDeriv[c] = new_value;
                if(fabs(new_value) > epsDeriv){
                    ok = true;
                }
            }

            // on ajoute la derivee a S:
            for(i=0 ; i<N ; i++){
                index = start[q] + dumMat(i, q);
                S[i] += clusterDeriv[index];
            }

        }
    }

    // Rprintf("other, nb iter=%i\n", iter);
    if(iter == iterMax){
        Rprintf("[Getting cluster deriv. other] Max iterations reached (%i)\n", iterMax);
    }

    return(S);
}

// Function to get the conditional sum of a matrix
[[cpp11::register]]
doubles_matrix<> cpp_tapply_sum(int Q, doubles_matrix<> x, integers dum){
    // Q: nber of classes
    // N: nber of observations
    // x: a matrix
    // dum: the N vector of clusters

    int N = x.nrow();
    int K = x.ncol();

    writable::doubles_matrix<> res(Q, K);
    int i, q, k;

    for(i=0 ; i<N ; i++){
        q = dum[i] - 1; // we take 1 off => different indexation in C

        for(k=0 ; k<K ; k++){
            res(q, k) += x(i, k);
        }
    }

    return(res);
}

// Function to get the conditional sum of a vector
[[cpp11::register]]
doubles cpp_tapply_vsum(int Q, doubles x, integers dum){
    // Q: nber of classes
    // x: a matrix
    // dum: the N vector of clusters

    int N = x.size();

    writable::doubles res(Q);
    int i, q;

    for(i=0 ; i<N ; i++){
        q = dum[i] - 1; // we take 1 off => different indexation in C
        res[q] += x[i];
    }

    return(res);
}

// similar a table but faster
[[cpp11::register]]
doubles cpp_table(int Q, integers dum){
    // Q: nber of classes
    // dum: the N vector of clusters

    int N = dum.size();

    writable::doubles res(Q);
    int i, q;

    for(i=0 ; i<N ; i++){
        q = dum[i] - 1; // we take 1 off => different indexation in C
        res[q]++;
    }

    return(res);
}

[[cpp11::register]]
double cpp_ssr_null(doubles y, doubles w = doubles(0)){
    // simple fun to compute the ssr of the null ols model
    // 2/3 times faster than pure r

    bool is_weight = w.size() > 1;

    int n = y.size();
    double denom = 0;

    // the "mean"
    double y_mean = 0;
    for(int i=0 ; i<n ; ++i){
        if(is_weight){
            y_mean += y[i] * w[i];
            denom += w[i];
        } else {
            y_mean += y[i];
        }
    }

    if(is_weight){
        y_mean = y_mean/denom;
    } else {
        y_mean = y_mean/n;
    }

    double res = 0, value = 0;
    for(int i=0 ; i<n ; ++i){
        value = y[i] - y_mean;
        if(is_weight){
            res += value * value * w[i];
        } else {
            res += value * value;
        }
    }

    return(res);
}

[[cpp11::register]]
double cpp_ssq(doubles x, doubles w = doubles(0)){
    // simple fun to compute the sum of the square of the elt of a vector
    // 30% faster than pure r (twice faster with weights)

    bool is_weight = w.size() > 1;

    int n = x.size();

    // the mean
    double res = 0;
    for(int i=0 ; i<n ; i++){
        if(is_weight){
            res += x[i] * x[i] * w[i];
        } else {
            res += x[i] * x[i];
        }
    }

    return res;
}

[[cpp11::register]]
bool cpp_isConstant(doubles x){
    // simple fun to see whether a variable is constant
    // it is unexpensive -- not the most useful function however:
    //		for 1e7 obs, you gain 50ms over var(x)==0, but it's still 50ms!

    int n = x.size();
    bool res = true;
    double value = x[0];
    for(int i=1 ; i<n ; i++){
        if(x[i] != value){
            res = false;
            break;
        }
    }

    return(res);
}

[[cpp11::register]]
bool cpp_any_na_null(SEXP x){
    // > twice faster than testing the two separately
    // x is a vector

    int n = Rf_length(x);
    double *px = REAL(x);

    for(int i=0 ; i<n ; ++i){
        double x_tmp = px[i];
        if(std::isnan(x_tmp) || x_tmp == 0){
            return true;
        }
    }

    return false;
}

[[cpp11::register]]
int cpp_constant_dum(int k, doubles x, integers dum, bool only_0 = false){
    // number of values of dum for which x is constant

    int n_obs = dum.size();

    double ref = x[0];
    int dum_current = dum[0];

    bool found_different = only_0 ? ref != 0 : false;
    int nb_constant = 0;

    for(int i=1 ; i<n_obs ; ++i){
        if(dum[i] != dum_current){
            // new guy
            dum_current = dum[i];
            if(found_different == false){
                ++nb_constant;
            }

            ref = x[i];
            found_different = only_0 ? ref != 0 : false;

        } else if(!found_different){
            if(x[i] != ref){
                found_different = true;
            }
        }

    }

    if(found_different == false){
        ++nb_constant;
    }

    return nb_constant;
}


//
// Lag related functions //
//

[[cpp11::register]]
list cpp_find_duplicates(integers id, integers time){
    // we check whether there are duplicated rows
    // if so, we provide information

    int n = id.size();
    int n_dup = 0;
    int obs_dup = 0;
    bool any_dup = false;

    /// finding the first duplicate value
    int i = 0;
    for(i=1 ; i<n ; ++i){
        if(time[i - 1] == time[i]){
            if(id[i - 1] == id[i]){
                any_dup = true;
                break;
            }
        }
    }

    // if dup: we find out the nber
    if(any_dup){
        obs_dup = i; // the 1 is implicitely added to make it r style
        int id_dup = id[i];
        int time_dup = time[i];
        n_dup = 2;
        while(++i<n && id_dup == id[i] && time_dup == time[i]) n_dup++;
    }

    writable::list res;
    res.push_back({"n_dup"_nm = n_dup});
    res.push_back({"obs_dup"_nm = obs_dup});

    return(res);
}

[[cpp11::register]]
int cpp_pgcd(integers x){
    // quick and dirty, but does not matter

    int n = x.size();

    if(n == 1){
        return(x[0]);
    }

    bool ok = false;
    int pgcd = x[0];

    // the min
    for(int i=1 ; i<n ; ++i){
        if(pgcd > x[i]){
            pgcd = x[i];
        }
    }

    // the denom
    while(!ok && pgcd > 1){
        ok = true;
        for(int i=0 ; i<n ; ++i){
            if(x[i] % pgcd != 0){
                pgcd--;
                ok = false;
                break;
            }
        }
    }


    return pgcd;
}

[[cpp11::register]]
integers cpp_lag_obs(integers id, integers time, int nlag){
    // in case of ties, we sum
    // must be two consecutive years
    // returns an observation nber of where to find the lagged obs
    // note that we forbid duplicate rows!!!!! => in case of duplicate we use a rule of thumb

    int nobs = id.size();
    int obs;
    int id_current;
    int time_current;

    //  OLD: res(nobs, NA_INTEGER)

    // NEW: pre-allocate vector with NA integer
    writable::integers res(nobs);
    for (int i = 0; i < nobs; i++) {
        res[i] = NA_INTEGER;
    }

    int i, j, diff_time;

    if(nlag > 0){
        i = 0;
        while(i < nobs){
            // R_CheckUserInterrupt(); // this is (too) costly
            id_current = id[i];
            time_current = time[i];
            obs = i + 1; // observation, R style
            j = i + 1;
            while(j < nobs){
                diff_time = time[j] - time_current;
                if(id[j] != id_current){
                    // we start someone else
                    i = j - 1; // minus 1 because after we indent i
                    break;
                } else if(diff_time > nlag){
                    // we are too far => stop
                    break;
                } else if(diff_time == 0){
                    // it's the same time/id => we skip
                    ++i;
                } else if(diff_time < nlag){
                    // not far enough
                } else if(diff_time == nlag){
                    // match!
                    res[j] = obs;
                }

                ++j;
            }
            ++i;
        }
    } else if(nlag < 0){
        /**************************************************************************
         * NOTA: I could have tweaked the previous if() to get rid of the condition
         *       but the code would have lost in clarity.
         *       For the lead: opposite to what is done before
         ***************************************************************************/
        int nlead = -nlag;
        i = nobs - 1;
        while(i >= 0){
            // R_CheckUserInterrupt(); // this is (too) costly
            id_current = id[i];
            time_current = time[i];
            obs = i + 1; // observation, R style
            j = i - 1;
            while(j >= 0){
                diff_time = time_current - time[j];
                if(id[j] != id_current){
                    // we start someone else
                    i = j + 1; // plus 1 because after we dedent i
                    break;
                } else if(diff_time > nlead){
                    // we are too far => stop
                    break;
                } else if(diff_time == 0){
                    // it's the same time/id => we skip
                    --i;
                } else if(diff_time < nlead){
                    // not far enough
                } else if(diff_time == nlead){
                    // match!
                    res[j] = obs;
                }

                --j;
            }
            --i;
        }
    } else {
        for(int i=0 ; i<nobs ; ++i){
            res[i] = i + 1;
        }
    }

    return(res);
}


[[cpp11::register]]
integers cpp_check_nested(SEXP fe_list, SEXP cluster_list, integers fe_sizes, int n){
    // Returns boolean vector of whether each FE is nested in the clusters

    int Q = Rf_length(fe_list);
    int G = Rf_length(cluster_list);

    // SEXP x0 = VECTOR_ELT(x, 0);

    writable::integers res(Q);

    for(int q=0 ; q<Q ; ++q){

        int *pfe = INTEGER(VECTOR_ELT(fe_list, q));

        for(int g=0 ; g<G ; ++g){
            vector<int> fe_clust(fe_sizes[q]);

            int *pclust =INTEGER(VECTOR_ELT(cluster_list, g));

            bool nested = true;
            int fe_value = 0;
            int clust_value = 0;
            for(int i=0 ; i<n ; ++i){
                fe_value = pfe[i] - 1;
                clust_value = fe_clust[fe_value];
                if(clust_value == 0){
                    fe_clust[fe_value] = pclust[i];
                } else if(clust_value != pclust[i]){
                    nested = false;
                    break;
                }
            }

            if(nested){
                res[q] = 1;
                break;
            }
        }
    }

    return res;
}


[[cpp11::register]]
doubles cpp_diag_XUtX(doubles_matrix<> X, doubles_matrix<> U){
    // computes the diagonal of X %*% U %*% t(X)

    int n = X.nrow();
    int K = X.ncol();

    writable::doubles res(n);

    for(int i=0 ; i<n ; ++i){

        double res_i = 0;
        for(int k=0 ; k<K ; ++k){

            double xk = 0;
            for(int k2=0 ; k2<K ; ++k2){
                xk += X(i, k2) * U(k, k2);
            }

            res_i += xk * X(i,k);
        }

        res[i] = res_i;
    }

    return res;
}


class simple_vec_double{
    simple_vec_double() = delete;
    double *px_double = nullptr;
    int *px_int = nullptr;
    int n;
    bool is_real;
public:
    simple_vec_double(SEXP x);
    double operator[](int);
};

simple_vec_double::simple_vec_double(SEXP x){

    n = Rf_length(x);

    if(TYPEOF(x) == REALSXP){
        px_double = REAL(x);
        is_real = true;

    } else if(TYPEOF(x) == INTSXP){
        px_int = INTEGER(x);
        is_real = false;

    } else {
        stop("Error: Wrong argument type in cpp_factor_matrix.");
    }
}

double simple_vec_double::operator[](int i){
    if(i >= n){
        return 1;
    } else if(is_real){
        return px_double[i];
    } else {
        return static_cast<double>(px_int[i]);
    }
}


[[cpp11::register]]
doubles_matrix<> cpp_factor_matrix(integers fact, logicals is_na_all, integers who_is_dropped, SEXP var, strings col_names){
    // fact: integer vector from 1 (!) to K, can contain NAs
    // Checking Na is cheap as opposed to populating the matrix, but having an argument avoids creating a new object

    int n = fact.size();
    int K = 0;

    // Finding out the TOTAL number of cols (before removal)
    for(int i=0 ; i<n ; ++i){
        if(!is_na_all[i] && K < fact[i]){
            K = fact[i];
        }
    }

    // Are there values to remove? If so, we create a mapping
    int n_drop = who_is_dropped[0] == -1 ? 0 : Rf_length(who_is_dropped);
    bool IS_REMOVAL = n_drop > 0;
    std::vector<int> mapping;

    if(IS_REMOVAL){
        mapping.resize(K);

        for(int k=0 ; k<K ; ++k){
            mapping[k] = k;
        }

        // In mapping: values equal to -1 mean dropped value

        int k_drop = 0;
        for(int k=0 ; k<K ; ++k){
            if(k_drop < n_drop && k + 1 == who_is_dropped[k_drop]){
                ++k_drop;
                mapping[k] = -1;
            } else {
                mapping[k] -= k_drop;
            }
        }

        K -= k_drop;
    }

    writable::doubles_matrix<> res(n, K);

    // The interacted var
    simple_vec_double my_var(var);

    // Filling the matrix
    for(int i=0 ; i<n ; ++i){
        if(is_na_all[i]){
            // we fill the row
            for(int k=0 ; k<K ; ++k){
                res(i, k) += NA_REAL;
            }

        } else if(IS_REMOVAL) {
            if(mapping[fact[i] - 1] != -1){
                res(i, mapping[fact[i] - 1]) = my_var[i];
            }

        } else {
            res(i, fact[i] - 1) = my_var[i];
        }
    }

    // TODO: dimnames
    // res.attr("dimnames") = list(NULL, col_names);

    return res;
}

[[cpp11::register]]
std::string cpp_add_commas(double x, int r = 1, bool whole = true){
    // a bit like (but not exactly equal to) format(x, nsmall = 1, big.mark = ",") but about 40-100 times faster
    // for whole numbers => no trailing digits
    // does not accept vectors, although super easy to expand to vectors

    std::string x_str = std::to_string(static_cast<int>(abs(x)));
    std::string res;

    if(x < 0){
        res.push_back('-');
        x = -x;
    }

    if(x < 1000){
        res.insert(res.size(), x_str);

    } else {
        int n = x_str.size();
        int e = n; // e: exponent

        while(e > 0){
            res.push_back(x_str[n - e]);
            --e;
            if(e > 1 && e % 3 == 0) {
                res.push_back(',');
            }
        }
    }

    double rest = x - floor(x);
    if((rest != 0 || !whole) && r > 0){
        // not a whole number

        res.push_back('.');

        if(r == 1){
            res.push_back(std::to_string(static_cast<int>(round(rest * 10)))[0]);

        } else {
            double rounded_rest = round(rest * pow(10, r)) / pow(10, r);
            std::string rest_str = std::to_string(rounded_rest);
            int nr = rest_str.size();

            for(int i=2 ; i < nr && i-1 <= r ; ++i){
                res.push_back(rest_str[i]);
            }
        }
    }

    return res;
}


[[cpp11::register]]
list cpp_find_never_always_treated(integers cohort, doubles period){
    // Note that both cohort and period are sorted according to cohort

    writable::integers always_treated;
    // cohort_ref => the guys that will be out of the interaction
    writable::integers cohort_ref;

    int n = cohort.size();

    bool is_pos = false;
    bool is_neg = false;
    bool is_ok = false;

    int j = cohort[0];
    if(period[0] < 0){
        is_neg = true;
    } else {
        is_pos = true;
    }

    for(int i=1 ; i<n ; ++i){

        if(j == cohort[i]){
            if(!is_ok){
                // condition avoids checking period
                // only worth for very (very) long vectors

                if(period[i] < 0){
                    is_neg = true;
                    is_ok = is_pos;
                } else {
                    is_pos = true;
                    is_ok = is_neg;
                }
            }
        } else {
            // we change IDs
            if(!is_ok){
                if(is_pos) always_treated.push_back(j);
                cohort_ref.push_back(j);
            }

            // re-init
            j = cohort[i];
            is_neg = false;
            is_pos = false;
            is_ok = false;

        }
    }

    // Last element
    if(!is_ok){
        if(is_pos) always_treated.push_back(j);
        cohort_ref.push_back(j);
    }

    writable::list res;
    res["always_treated"] = always_treated;
    res.push_back({"ref"_nm = cohort_ref});

    return res;
}


[[cpp11::register]]
integers cpp_get_first_item(integers x, int n_items){
    // observation id of the first occurrence
    // x ranges from 1 to n_items
    // we return indexes R style

    int n = x.size();
    writable::integers res(n_items);

    for(int i=0 ; i<n ; ++i){
        if(res[x[i] - 1] == 0){
            res[x[i] - 1] = i + 1;
        }
    }

    return res;
}


[[cpp11::register]]
integers cpp_combine_clusters(SEXP cluster_list, integers index){
    // cluster: list of integer vectors, each ranging from 1 to the number of cases
    // index: result of order() on the clusters

    if(TYPEOF(cluster_list) != VECSXP){
        stop("Internal error: Only lists are accepted!");
    }

    int Q = Rf_length(cluster_list);

    int n = index.size();
    writable::integers res(n);

    // Loading the data
    vector<int*> pcluster(Q);
    for(int q=0 ; q<Q ; ++q){
        SEXP cluster_q = VECTOR_ELT(cluster_list, q);

        pcluster[q] = INTEGER(cluster_q);
    }

    // the observation ID
    int obs = index[0] - 1;

    // vector holding the current value
    vector<int> current_value(Q);

    // initialization
    int counter = 1;
    res[obs] = counter;
    for(int q=0 ; q<Q ; ++q){
        current_value[q] = pcluster[q][obs];
    }

    // we loop on the vector and flag values that are different
    int q = 0;
    for(int i=1 ; i<n ; ++i){
        obs = index[i] - 1;

        for(q=0 ; q<Q ; ++q){
            if(pcluster[q][obs] != current_value[q]){
                break;
            }
        }

        // if the condition holds => means the values are different
        if(q < Q){
            ++counter;
            // we save the new values
            for(; q<Q ; ++q){
                current_value[q] = pcluster[q][obs];
            }
        }

        res[obs] = counter;
    }

    return res;
}




[[cpp11::register]]
list cpp_cut(doubles x_sorted, doubles cut_points, integers is_included){
    // x_sorted: no NA, sorted
    // cut_points: bounds
    // is_included: for each bound, if it is included or not
    //

    int N = x_sorted.size();
    int n_cuts = cut_points.size();

    bool is_int = true;
    for(int i=0 ; i<N ; ++i){
        if(fabs(x_sorted[i] - round(x_sorted[i])) > 0.00000000001){
            is_int = false;
            break;
        }
    }

    writable::integers x_int(N);
    for (int i=0; i < N; i++) {
        x_int[i] = n_cuts + 1;
    }

    writable::integers isnt_empty(n_cuts + 1);
    writable::doubles value_min(n_cuts + 1);
    writable::doubles value_max(n_cuts + 1);

    int index = 0;
    bool first = true;
    bool include = is_included[0];
    double cutoff = cut_points[0];
    int i = 0;
    while(i < N){

        if(include ? x_sorted[i] <= cutoff : x_sorted[i] < cutoff){

            if(first){
                isnt_empty[index] = true;
                value_min[index] = x_sorted[i];
                first = false;
            }

            x_int[i] = index + 1;
            ++i;

        } else {
            // we increment index
            // bins can be empty

            if(isnt_empty[index] && i > 0){
                value_max[index] = x_sorted[i - 1];
            }

            ++index;

            if(index == n_cuts){
                // last bin: we've done it at initialization
                isnt_empty[index] = true;
                value_min[index] = x_sorted[i];
                value_max[index] = x_sorted[N - 1];
                break;
            }

            include = is_included[index];
            cutoff = cut_points[index];
            first = true;
        }
    }

    if(index != n_cuts){
        value_max[index] = x_sorted[N - 1];
    }

    writable::list res;
    res["x_int"] = x_int;
    res.push_back({"isnt_empty"_nm = isnt_empty});
    res.push_back({"value_min"_nm = value_min});
    res.push_back({"value_max"_nm = value_max});
    res.push_back({"is_int"_nm = is_int});

    return res;
}

[[cpp11::register]]
bool cpp_is_int(SEXP x){

    if(TYPEOF(x) == INTSXP){
        return true;
    }

    if(TYPEOF(x) != REALSXP){
        return false;
    }

    int N = Rf_length(x);
    double *px = REAL(x);

    bool is_int = true;
    for(int i=0 ; i<N ; ++i){
        if(fabs(px[i] - round(px[i])) > 0.00000000001){
            is_int = false;
            break;
        }
    }

    return is_int;
}

[[cpp11::register]]
double cpp_hash_string(std::string x){
    // simple function hashing a string
    // used to identify tables in etable

    double res = std::hash<std::string>{}(x);

    return res;
}
