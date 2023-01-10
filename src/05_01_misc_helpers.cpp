#include "05_0_misc.h"

[[cpp11::register]]
std::string cpp_add_commas(double x, int r = 1, bool whole = true)
{
    // a bit like (but not exactly equal to) format(x, nsmall = 1, big.mark = ",") but about 40-100 times faster
    // for whole numbers => no trailing digits
    // does not accept vectors, although super easy to expand to vectors

    std::string x_str = std::to_string(static_cast<int>(abs(x)));
    std::string res;

    if (x < 0)
    {
        res.push_back('-');
        x = -x;
    }

    if (x < 1000)
    {
        res.insert(res.size(), x_str);
    }
    else
    {
        int n = x_str.size();
        int e = n; // e: exponent

        while (e > 0)
        {
            res.push_back(x_str[n - e]);
            --e;
            if (e > 1 && e % 3 == 0)
            {
                res.push_back(',');
            }
        }
    }

    double rest = x - floor(x);
    if ((rest != 0 || !whole) && r > 0)
    {
        // not a whole number

        res.push_back('.');

        if (r == 1)
        {
            res.push_back(std::to_string(static_cast<int>(round(rest * 10)))[0]);
        }
        else
        {
            double rounded_rest = round(rest * pow(10, r)) / pow(10, r);
            std::string rest_str = std::to_string(rounded_rest);
            int nr = rest_str.size();

            for (int i = 2; i < nr && i - 1 <= r; ++i)
            {
                res.push_back(rest_str[i]);
            }
        }
    }

    return res;
}

[[cpp11::register]] list cpp_find_never_always_treated(integers cohort, doubles period)
{
    // Note that both cohort and period are sorted according to cohort

    writable::integers always_treated;
    // cohort_ref => the guys that will be out of the interaction
    writable::integers cohort_ref;

    int n = cohort.size();

    bool is_pos = false;
    bool is_neg = false;
    bool is_ok = false;

    int j = cohort[0];
    if (period[0] < 0)
    {
        is_neg = true;
    }
    else
    {
        is_pos = true;
    }

    for (int i = 1; i < n; ++i)
    {

        if (j == cohort[i])
        {
            if (!is_ok)
            {
                // condition avoids checking period
                // only worth for very (very) long vectors

                if (period[i] < 0)
                {
                    is_neg = true;
                    is_ok = is_pos;
                }
                else
                {
                    is_pos = true;
                    is_ok = is_neg;
                }
            }
        }
        else
        {
            // we change IDs
            if (!is_ok)
            {
                if (is_pos)
                    always_treated.push_back(j);
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
    if (!is_ok)
    {
        if (is_pos)
            always_treated.push_back(j);
        cohort_ref.push_back(j);
    }

    writable::list res;
    res.push_back({"always_treated"_nm = always_treated});
    res.push_back({"ref"_nm = cohort_ref});

    return res;
}

[[cpp11::register]] integers cpp_get_first_item(integers x, int n_items)
{
    // observation id of the first occurrence
    // x ranges from 1 to n_items
    // we return indexes R style

    int n = x.size();
    writable::integers res(n_items);

    for (int i = 0; i < n; ++i)
    {
        if (res[x[i] - 1] == 0)
        {
            res[x[i] - 1] = i + 1;
        }
    }

    return res;
}

[[cpp11::register]] integers cpp_combine_clusters(SEXP cluster_list, integers index)
{
    // cluster: list of integer vectors, each ranging from 1 to the number of cases
    // index: result of order() on the clusters

    if (TYPEOF(cluster_list) != VECSXP)
    {
        stop("Internal error: Only lists are accepted!");
    }

    int Q = Rf_length(cluster_list);

    int n = index.size();
    writable::integers res(n);

    // Loading the data
    vector<int *> pcluster(Q);
    for (int q = 0; q < Q; ++q)
    {
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
    for (int q = 0; q < Q; ++q)
    {
        current_value[q] = pcluster[q][obs];
    }

    // we loop on the vector and flag values that are different
    int q = 0;
    for (int i = 1; i < n; ++i)
    {
        obs = index[i] - 1;

        for (q = 0; q < Q; ++q)
        {
            if (pcluster[q][obs] != current_value[q])
            {
                break;
            }
        }

        // if the condition holds => means the values are different
        if (q < Q)
        {
            ++counter;
            // we save the new values
            for (; q < Q; ++q)
            {
                current_value[q] = pcluster[q][obs];
            }
        }

        res[obs] = counter;
    }

    return res;
}

[[cpp11::register]] list cpp_cut(doubles x_sorted, doubles cut_points, integers is_included)
{
    // x_sorted: no NA, sorted
    // cut_points: bounds
    // is_included: for each bound, if it is included or not
    //

    int N = x_sorted.size();
    int n_cuts = cut_points.size();

    bool is_int = true;
    for (int i = 0; i < N; ++i)
    {
        if (fabs(x_sorted[i] - round(x_sorted[i])) > 0.00000000001)
        {
            is_int = false;
            break;
        }
    }

    writable::integers x_int(N);
    for (int i = 0; i < N; i++)
    {
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
    while (i < N)
    {

        if (include ? x_sorted[i] <= cutoff : x_sorted[i] < cutoff)
        {

            if (first)
            {
                isnt_empty[index] = true;
                value_min[index] = x_sorted[i];
                first = false;
            }

            x_int[i] = index + 1;
            ++i;
        }
        else
        {
            // we increment index
            // bins can be empty

            if (isnt_empty[index] && i > 0)
            {
                value_max[index] = x_sorted[i - 1];
            }

            ++index;

            if (index == n_cuts)
            {
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

    if (index != n_cuts)
    {
        value_max[index] = x_sorted[N - 1];
    }

    writable::list res;
    res.push_back({"x_int"_nm = x_int});
    res.push_back({"isnt_empty"_nm = isnt_empty});
    res.push_back({"value_min"_nm = value_min});
    res.push_back({"value_max"_nm = value_max});
    res.push_back({"is_int"_nm = is_int});

    return res;
}

[[cpp11::register]] bool cpp_is_int(SEXP x)
{

    if (TYPEOF(x) == INTSXP)
    {
        return true;
    }

    if (TYPEOF(x) != REALSXP)
    {
        return false;
    }

    int N = Rf_length(x);
    double *px = REAL(x);

    bool is_int = true;
    for (int i = 0; i < N; ++i)
    {
        if (fabs(px[i] - round(px[i])) > 0.00000000001)
        {
            is_int = false;
            break;
        }
    }

    return is_int;
}

[[cpp11::register]] double cpp_hash_string(std::string x)
{
    // simple function hashing a string
    // used to identify tables in etable

    double res = std::hash<std::string>{}(x);

    return res;
}

[[cpp11::register]] bool cpp_isConstant(doubles x)
{
    // simple fun to see whether a variable is constant
    // it is unexpensive -- not the most useful function however:
    //		for 1e7 obs, you gain 50ms over var(x)==0, but it's still 50ms!

    int n = x.size();
    bool res = true;
    double value = x[0];
    for (int i = 1; i < n; i++)
    {
        if (x[i] != value)
        {
            res = false;
            break;
        }
    }

    return (res);
}

[[cpp11::register]] bool cpp_any_na_null(SEXP x)
{
    // > twice faster than testing the two separately
    // x is a vector

    int n = Rf_length(x);
    double *px = REAL(x);

    for (int i = 0; i < n; ++i)
    {
        double x_tmp = px[i];
        if (std::isnan(x_tmp) || x_tmp == 0)
        {
            return true;
        }
    }

    return false;
}

[[cpp11::register]] list cpp_find_duplicates(integers id, integers time)
{
    // we check whether there are duplicated rows
    // if so, we provide information

    int n = id.size();
    int n_dup = 0;
    int obs_dup = 0;
    bool any_dup = false;

    /// finding the first duplicate value
    int i = 0;
    for (i = 1; i < n; ++i)
    {
        if (time[i - 1] == time[i])
        {
            if (id[i - 1] == id[i])
            {
                any_dup = true;
                break;
            }
        }
    }

    // if dup: we find out the nber
    if (any_dup)
    {
        obs_dup = i; // the 1 is implicitely added to make it r style
        int id_dup = id[i];
        int time_dup = time[i];
        n_dup = 2;
        while (++i < n && id_dup == id[i] && time_dup == time[i])
            n_dup++;
    }

    writable::list res;
    res.push_back({"n_dup"_nm = n_dup});
    res.push_back({"obs_dup"_nm = obs_dup});

    return (res);
}

[[cpp11::register]] integers cpp_check_nested(SEXP fe_list, SEXP cluster_list, integers fe_sizes, int n)
{
    // Returns boolean vector of whether each FE is nested in the clusters

    int Q = Rf_length(fe_list);
    int G = Rf_length(cluster_list);

    // SEXP x0 = VECTOR_ELT(x, 0);

    writable::integers res(Q);

    for (int q = 0; q < Q; ++q)
    {

        int *pfe = INTEGER(VECTOR_ELT(fe_list, q));

        for (int g = 0; g < G; ++g)
        {
            vector<int> fe_clust(fe_sizes[q]);

            int *pclust = INTEGER(VECTOR_ELT(cluster_list, g));

            bool nested = true;
            int fe_value = 0;
            int clust_value = 0;
            for (int i = 0; i < n; ++i)
            {
                fe_value = pfe[i] - 1;
                clust_value = fe_clust[fe_value];
                if (clust_value == 0)
                {
                    fe_clust[fe_value] = pclust[i];
                }
                else if (clust_value != pclust[i])
                {
                    nested = false;
                    break;
                }
            }

            if (nested)
            {
                res[q] = 1;
                break;
            }
        }
    }

    return res;
}

[[cpp11::register]] int cpp_pgcd(integers x)
{
    // quick and dirty, but does not matter

    int n = x.size();

    if (n == 1)
    {
        return (x[0]);
    }

    bool ok = false;
    int pgcd = x[0];

    // the min
    for (int i = 1; i < n; ++i)
    {
        if (pgcd > x[i])
        {
            pgcd = x[i];
        }
    }

    // the denom
    while (!ok && pgcd > 1)
    {
        ok = true;
        for (int i = 0; i < n; ++i)
        {
            if (x[i] % pgcd != 0)
            {
                pgcd--;
                ok = false;
                break;
            }
        }
    }

    return pgcd;
}

//
// Getting the cluster coefficients
//

[[cpp11::register]]
list cpp_get_fe_gnl(int Q, int N, writable::doubles sumFE, writable::integers_matrix<> dumMat, writable::integers cluster_sizes, writable::integers obsCluster)
{
    // This function returns a list of the cluster coefficients for each cluster
    // dumMat: the matrix of cluster ID for each observation, with cpp index style
    // Q, N: nber of clusters / obs
    // cluster_sizes: vector of the number of cases per cluster
    // obsCluster: the integer vector that is equal to order(dum[[g]])
    // RETURN:
    // a list for each cluster of the cluster coefficient value
    // + the last element is the number of clusters that have been set as references (nb_ref)
    // => much optimized version May 2019

    int iter = 0, iterMax = 10000;
    int iter_loop = 0, iterMax_loop = 10000;

    // Creation of the indices to put all the cluster values into a single vector
    int nb_coef = 0;
    int nb_ref[Q]; // nb_ref takes the nb of elements set as ref

    for (int q = 0; q < Q; q++)
    {
        // the total number of clusters (eg if c1: man/woman and c2: 10 countries: total of 12 cases)
        nb_coef += cluster_sizes[q];
    }

    double cluster_values[nb_coef];
    int cluster_visited[nb_coef]; // whether a value has been already assigned

    // index of the cluster
    std::vector<int*> pindex_cluster(Q);
    std::vector<int> index_cluster(nb_coef);

    for (int i = 0; i < nb_coef; i++)
    {
        index_cluster[i] += i;
    }

    // TODO: .data() comes from Rcpp, what does it exactly do?
    pindex_cluster[0] = index_cluster.data();
    for (int q = 1; q < Q; q++)
    {
        pindex_cluster[q] = pindex_cluster[q - 1] + cluster_sizes[q - 1];
    }

    // Now we create the vector of observations for each cluster
    // we need a starting and an end vector as well
    int start_cluster[nb_coef];
    int end_cluster[nb_coef];

    int index;
    int k;

    for (int q = 0; q < Q; q++)
    {
        // table cluster: nber of elements for each cluster class
        writable::integers tableCluster;
        tableCluster.push_back(cluster_sizes[q]);
        for (int i = 0; i < N; i++)
        {
            k = dumMat(i, q);
            tableCluster[k] += 1; // the number of obs per case
        }

        // now creation of the start/end vectors
        for (int k = 0; k < cluster_sizes[q]; k++)
        {
            int *pindex = pindex_cluster[q];
            index = pindex[k];
            // index = start(q) + k;
            if (k == 0) {
                start_cluster[index] = 0;
                end_cluster[index] = tableCluster[k];
            } else {
                start_cluster[index] = end_cluster[index - 1];
                end_cluster[index] = end_cluster[index - 1] + tableCluster[k];
            }
        }
    }

    // matrix of the clusters that have been computed
    writable::integers_matrix<> mat_done(N, Q);
    writable::integers rowsums(N);

    // vector of elements to loop over
    writable::integers id2do(N);
    writable::integers id2do_next(N);
    int nb2do = N, nb2do_next = N;
    for (int i = 0; i < nb2do; i++)
    {
        id2do[i] = i;
        id2do_next[i] = i;
    }

    // Other indices and variables
    int qui_max, obs;
    int rs, rs_max; // rs: row sum
    int id_cluster;
    double other_value;
    bool first;

    //
    // THE MAIN LOOP
    //

    while (iter < iterMax)
    {
        iter++;

        //
        // Finding the row where to put the 0s
        //

        if (iter == 1)
        {
            // 1st iter, we select the first element
            qui_max = 0;
        }
        else
        {
            // we find the row that has the maximum of items done

            qui_max = 0;
            rs_max = 0;
            for (int i = 0; i < nb2do; i++)
            {
                obs = id2do[i];

                rs = rowsums[obs];

                if (rs == Q - 2)
                {
                    // if rs == Q-2 => its the maximum possible, no need to go further
                    qui_max = obs;
                    break;
                }
                else if (rs < Q && rs > rs_max)
                {
                    // this is to handle complicated cases with more than 3+ clusters
                    qui_max = obs;
                    rs_max = rs;
                }
                else if (qui_max == 0 && rs == 0)
                {
                    // used to initialize qui_max
                    qui_max = obs;
                }
            }
        }

        // Rprintf("Obs selected: %i\n", qui_max);

        //
        // Putting the 0s, ie setting the references
        //

        // the first element is spared
        first = true;
        for (int q = 0; q < Q; q++)
        {
            if (mat_done(qui_max, q) == 0)
            {
                if (first)
                {
                    // we spare the first element
                    first = false;
                }
                else
                {
                    // we set the cluster to 0
                    // 1) we find the cluster
                    id_cluster = dumMat(qui_max, q);
                    // Rprintf("Cluster: %i\n", id_cluster + 1);
                    // 2) we get the index of the cluster vector
                    // TODO: no viable conversion from integers to int* here
                    int *pindex = pindex_cluster[q];
                    index = pindex[id_cluster];
                    // 3) we set the cluster value to 0
                    cluster_values[index] = 0;
                    // 4) we update the mat_done matrix for the elements of this cluster
                    // TODO: int/integers conversion
                    for (int i = start_cluster[index]; i < end_cluster[index]; i++)
                    {
                        // TODO: why does this subset a vector in a matrix way??
                        // original https://github.com/pachadotdev/fixest2/blob/master/src/misc_funs.cpp#L382
                        obs = obsCluster[i, q];
                        mat_done(obs, q) = 1;
                        rowsums[obs]++;
                    }
                    // 5) => we save the information on which cluster was set as a reference
                    nb_ref[q]++;
                }
            }
        }

        //
        // LOOP OF ALL OTHER UPDATES (CRITICAL)
        //

        iter_loop = 0;
        while (iter_loop < iterMax_loop)
        {
            iter_loop++;

            // Rprintf("nb2do_next: %i -- nb2do: %i\n", nb2do_next, nb2do);

            R_CheckUserInterrupt();

            //
            // Selection of indexes (new way) to be updated
            //

            // initialisation of the observations to cover (before first loop the two are identical)
            if (iter_loop != 1)
            {
                nb2do = nb2do_next;
                for (int i = 0; i < nb2do; i++)
                {
                    id2do[i] += id2do_next[i];
                }
            }

            nb2do_next = 0;

            for (int i = 0; i < nb2do; i++)
            {
                // we compute the rowsum of the obs that still needs to be done
                obs = id2do[i];

                rs = rowsums[obs];

                if (rs < Q - 1)
                {
                    // you need to check it next time!
                    id2do_next[nb2do_next] = obs;
                    nb2do_next++;
                }
                else if (rs == Q - 1)
                {
                    // means: needs to be updated
                    int q;
                    for (q = 0; mat_done(obs, q) != 0; q++)
                    {
                        // selection of the q that is equal to 0
                    }

                    int *pindex = pindex_cluster[q];
                    int index_select = pindex[dumMat(obs, q)];

                    // Finding the value of the cluster coefficient
                    other_value = 0.0;
                    // Computing the sum of the other cluster values
                    // and finding the cluster to be updated (q_select)
                    for (int l = 0; l < Q; l++)
                    {
                        // we can loop over all q because cluster_values is initialized to 0
                        // TODO: go back and convert pindex_cluster to matrix?
                        index = pindex_cluster[l][dumMat(obs, l)];
                        // TODO: double vs cpp11::doubles?
                        other_value += cluster_values[index];
                    }

                    // the index to update
                    // TODO: no viable overloaded +=
                    cluster_values[index_select] += sumFE[obs] - other_value;

                    // Update of the mat_done
                    for (int j = start_cluster[index_select]; j < end_cluster[index_select]; j++)
                    {
                        // TODO: comma operator here is meaningless
                        // original: https://github.com/pachadotdev/fixest2/blob/master/src/misc_funs.cpp#L453
                        obs = obsCluster[j, q];
                        mat_done(obs, q) = 1;
                        rowsums[obs]++;
                    }
                }
            }

            if (nb2do_next == nb2do)
                break;
        }

        // out of this loop:
        //	- nb2do == nb2do_next
        // - id2do == id2do_next

        // Check that everything is all right
        if (iter_loop == iterMax_loop)
        {
            Rprintf("Problem getting FE, maximum iterations reached (2nd order loop).");
        }

        // if no more obs to be covered
        if (nb2do_next == 0)
            break;
    }

    if (iter == iterMax)
    {
        Rprintf("Problem getting FE, maximum iterations reached (1st order loop).");
    }

    // final formatting and save
    writable::list res(Q + 1);

    int *pindex;
    for (int q = 0; q < Q; q++)
    {
        writable::doubles quoi(cluster_sizes[q]);
        // TODO: incompatible integers vs int*
        pindex = pindex_cluster[q];
        for (k = 0; k < cluster_sizes[q]; k++)
        {
            index = pindex[k];
            // index = start(q) + k;
            // TODO: overload
            quoi[k] = cluster_values[index];
        }
        res[q] = quoi;
    }

    // TODO: viable overload
    res[Q] = nb_ref;

    return res;
}
