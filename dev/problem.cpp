#include <cpp11.hpp>
#include <numeric>
#include <vector> // should I just use regular arrays?
#include <iostream>

using namespace cpp11;
using namespace std;

[[cpp11::register]] list cpp_get_fe_gnl_(int Q,
                                        int N,
                                        writable::doubles sumFE,
                                        writable::integers_matrix<> dumMat,
                                        writable::integers cluster_sizes,
                                        writable::integers_matrix<> obsCluster)
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
    writable::integers nb_ref(Q); // nb_ref takes the nb of elements set as ref

    for (int q = 0; q < Q; q++)
    {
        // the total number of clusters (eg if c1: man/woman and c2: 10 countries: total of 12 cases)
        nb_coef = cluster_sizes[q];
    }

    writable::doubles cluster_values(nb_coef);
    writable::integers cluster_visited(nb_coef);

    // index of the cluster
    std::vector<int*> pindex_cluster(Q);
    std::vector<int> index_cluster(nb_coef);

    for (int i = 0; i < nb_coef; i++)
    {
        index_cluster[i] = i;
    }

    pindex_cluster[0] = index_cluster.data();
    for (int q = 1; q < Q; q++)
    {
        pindex_cluster[q] = pindex_cluster[q - 1] + cluster_sizes[q - 1];
    }

    // Now we create the vector of observations for each cluster
    // we need a starting and an end vector as well

    // writable::integers start_cluster(nb_coef);
    // writable::integers end_cluster(nb_coef);
    // => too much blah blah with operator=

    // maybe define start_cluster as an array instead
    // int start_cluster[nb_coef];
    // int end_cluster[nb_coef];
    // => warning: ISO C++ forbids variable length array ‘start_cluster’ [-Wvla]

    // alternative?
    std::vector<int> start_cluster(nb_coef);
    std::vector<int> end_cluster(nb_coef);

    int index;
    int k;

    for (int q = 0; q < Q; q++)
    {
        // table cluster: nber of elements for each cluster class
        writable::integers tableCluster(cluster_sizes.at(q));

        for (int i = 0; i < N; i++)
        {
            k = dumMat(i, q);
            tableCluster[k] = 1; // the number of obs per case
        }

        // now creation of the start/end vectors
        for (int k = 0; k < cluster_sizes[q]; k++)
        {
            int *pindex = pindex_cluster[q];
            index = pindex[k];

            // index = start(q) + k;
            if (k == 0)
            {
                start_cluster[index] = 0;
                end_cluster[index] = tableCluster[k];
            }
            else
            {
                start_cluster[index] = end_cluster[index - 1];
                end_cluster[index] = end_cluster[index - 1] + tableCluster[k];
            }
        }
    }
    
    // matrix of the clusters that have been computed
    writable::integers_matrix<> mat_done(N, Q);
    writable::integers rowsums(N);

    // fill with 0s
    for (int i = 0; i < N; i++)
    {
        rowsums[i] = 0;
        for (int j = 0; j < Q; j++)
        {
            mat_done(i, j) = 0;
        }
    }

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

        // Finding the row where to put the 0s

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

        // Putting the 0s, ie setting the references

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
                    // 2) we get the index of the cluster vector
                    int *pindex = pindex_cluster[q];
                    index = pindex[id_cluster];
                    // 3) we set the cluster value to 0
                    cluster_values[index] = 0;
                    // 4) we update the mat_done matrix for the elements of this cluster
                    for (int i = start_cluster[index]; i < end_cluster[index]; i++)
                    {
                        obs = obsCluster(i, q);
                        mat_done(obs, q) = 1;
                        rowsums[obs]++;
                    }
                    // 5) => we save the information on which cluster was set as a reference
                    nb_ref[q]++;
                }
            }
        }

        // LOOP OF ALL OTHER UPDATES (CRITICAL)

        iter_loop = 0;
        while (iter_loop < iterMax_loop)
        {
            iter_loop++;

            R_CheckUserInterrupt();

            // Selection of indexes (new way) to be updated

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
                        index = pindex_cluster[l][dumMat(obs, l)];
                        other_value += cluster_values[index];
                    }

                    // the index to update
                    cluster_values[index_select] += sumFE[obs] - other_value;

                    // Update of the mat_done
                    for (int j = start_cluster[index_select]; j < end_cluster[index_select]; j++)
                    {
                        obs = obsCluster(j, q);
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
        if (nb2do_next == 0) break;
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
        // writable::doubles quoi(cluster_sizes[q]);
        // same blah blah with operator=
        // L865 will say error: ambiguous overload for ‘operator=’
        // (operand types are ‘cpp11::writable::r_vector<double>::proxy’ and 
        // ‘cpp11::writable::r_vector<double>::proxy’)

        writable::doubles quoi(cluster_sizes[q]);

        pindex = pindex_cluster[q];
        for (k = 0; k < cluster_sizes[q]; k++)
        {
            index = pindex[k];
            // error: ambiguous overload for ‘operator=’ (operand types are ‘cpp11::writable::r_vector<double>::proxy’
            // quoi[k] = cluster_values[index];
            // why does it have to be like this?
            // we previously had writable::doubles cluster_values(nb_coef);
            // quoi[k] = static_cast<double>(cluster_values[index]);
            // quoi[k] = cluster_values[index];
            quoi.operator[](k).operator=(cluster_values.operator[](index));
        }
        res[q] = quoi;
    }

    res[Q] = nb_ref;

    return res;
}
