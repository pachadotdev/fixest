#include "06_0_quf.h"

inline unsigned long long float_to_ull(void *u, int i)
{
    unsigned long long *pu_ll = reinterpret_cast<unsigned long long *>(u);
    unsigned long long u_ull = pu_ll[i];
    unsigned long long mask = -(u_ull >> 63) | 0x8000000000000000;
    return (u_ull ^ mask);
}

double ull_to_float(unsigned long long u_ull)
{
    unsigned long long mask = ((u_ull >> 63) - 1) | 0x8000000000000000;
    unsigned long long res_ull = u_ull ^ mask;
    unsigned long long *pres_ull = &res_ull;
    double *pres = reinterpret_cast<double *>(pres_ull);
    double res = *pres;
    return res;
}

void quf_double(int n, int *x_uf, void *px, vector<double> &x_unik, bool is_string = false)
{

    // x_uf: x unclassed factor
    // px: pointer to x vector (R vector) -- READ ONLY!!!
    // x_unik: empty vector
    // px: either double or ULL (in case of strings)

    double *px_dble = (double *)px;
    unsigned long long *px_ull = (unsigned long long *)px;

    // variables
    unsigned long long x_ull_current = 0;
    // char: 2 octets: 0-255
    // one double is made of 8 char
    unsigned char x_char = 0;
    vector<unsigned long long> x_ulong(is_string ? 1 : n), x_tmp(n);

    // in case the data is string, we use px as the x_ulong
    unsigned long long *px_ulong = is_string ? px_ull : x_ulong.data();

    // radix array
    int radix_table[8][256] = {{0}};

    // 1) Counting
    for (int i = 0; i < n; ++i)
    {

        if (is_string)
        {
            x_ull_current = px_ull[i];
        }
        else
        {
            // we change x double into ulong after flipping to keep order
            x_ull_current = float_to_ull(px_dble, i);
            x_ulong[i] = x_ull_current;
        }

        for (int b = 0; b < 8; ++b)
        {
            x_char = ((unsigned char *)&x_ull_current)[b];
            ++radix_table[b][x_char];
        }
    }

    // 1') skipping
    vector<bool> skip_flag(8);
    for (int b = 0; b < 8; ++b)
    {
        x_char = ((unsigned char *)&x_ull_current)[b];
        skip_flag[b] = radix_table[b][x_char] == n;
    }

    // 2) Cumulating
    for (int b = 0; b < 8; ++b)
    {
        for (int d = 1; d < 256; ++d)
        {
            radix_table[b][d] += radix_table[b][d - 1];
        }
    }

    // 3) Sorting
    unsigned long long *x_read = px_ulong;
    unsigned long long *x_write = x_tmp.data();
    unsigned long long *p_tmp; // used for swapping

    vector<int> x_order(n), x_unclass(n);
    int *o_read = x_unclass.data();
    int *o_write = x_order.data();
    int *o_tmp;

    int k = 0;
    for (auto &p : x_unclass)
        p = k++;

    for (int b = 0; b < 8; ++b)
    {
        if (skip_flag[b] == false)
        {
            int index;
            for (int i = n - 1; i >= 0; --i)
            {
                x_ull_current = x_read[i];
                x_char = ((unsigned char *)&x_ull_current)[b];
                index = --radix_table[b][x_char];
                x_write[index] = x_read[i];
                o_write[index] = o_read[i];
            }

            p_tmp = x_read;
            x_read = x_write;
            x_write = p_tmp;

            o_tmp = o_read;
            o_read = o_write;
            o_write = o_tmp;
        }
    }

    if (o_read != x_order.data())
    {
        // we copy the result
        memcpy(x_order.data(), o_read, sizeof(int) * n);
    }

    // We unclass, starting at 1
    k = 1;
    x_unclass[0] = k;

    unsigned long long xi, xim1 = x_read[0];
    x_unik.push_back(is_string ? static_cast<double>(x_order[0] + 1) : ull_to_float(xim1));

    for (int i = 1; i < n; ++i)
    {
        xi = x_read[i];
        if (xi != xim1)
        {
            ++k;
            x_unik.push_back(is_string ? static_cast<double>(x_order[i] + 1) : ull_to_float(xi));
        }
        x_unclass[i] = k;
        xim1 = xi;
    }

    // We put into the right order
    for (int i = 0; i < n; ++i)
    {
        x_uf[x_order[i]] = x_unclass[i];
    }
}

void quf_int_gnl(int n, int *x_uf, void *px, vector<double> &x_unik, int x_min, bool is_double)
{
    // we can sort a range up to 2**31 => ie not the full int range
    // for ranges > 2**31 => as double
    // px: pointer to the values of x -- R vector READ ONLY!!!
    // px: either a R vector of type integer // either a R vector of type double
    // px: we know it thanks to the value of 'is_double'

    int *px_int = (int *)px;
    double *px_dble = (double *)px;

    // variables
    int x_uint_current = 0;
    vector<int> x_uint(n), x_tmp(n);

    // radix array
    int radix_table[4][256] = {{0}};

    // 1) Counting
    for (int i = 0; i < n; ++i)
    {
        // x_uint_current = px[i] - x_min;
        x_uint_current = is_double ? static_cast<int>(px_dble[i] - x_min) : px_int[i] - x_min;
        x_uint[i] = x_uint_current;

        for (int b = 0; b < 4; ++b)
        {
            ++radix_table[b][(x_uint_current >> 8 * b) & 0xFF];
        }
    }

    // 1') skipping
    vector<bool> skip_flag(4);
    for (int b = 0; b < 4; ++b)
    {
        skip_flag[b] = radix_table[b][(x_uint_current >> 8 * b) & 0xFF] == n;
    }

    // 2) Cumulating
    for (int b = 0; b < 4; ++b)
    {
        for (int d = 1; d < 256; ++d)
        {
            radix_table[b][d] += radix_table[b][d - 1];
        }
    }

    // 3) Sorting
    int *x_read = x_uint.data();
    int *x_write = x_tmp.data();
    int *p_tmp; // used for swapping

    vector<int> x_order(n), x_unclass(n);
    int *o_read = x_unclass.data();
    int *o_write = x_order.data();
    int *o_tmp;

    int k = 0;
    for (auto &p : x_unclass)
        p = k++;

    for (int b = 0; b < 4; ++b)
    {
        if (skip_flag[b] == false)
        {
            // Rprintf("bin: %i\n", b);
            int index;
            for (int i = n - 1; i >= 0; --i)
            {
                index = --radix_table[b][(x_read[i] >> 8 * b) & 0xFF];
                x_write[index] = x_read[i];
                o_write[index] = o_read[i];
            }

            p_tmp = x_read;
            x_read = x_write;
            x_write = p_tmp;

            o_tmp = o_read;
            o_read = o_write;
            o_write = o_tmp;
        }
    }

    if (o_read != x_order.data())
    {
        // we copy the result, if needed
        memcpy(x_order.data(), o_read, sizeof(int) * n);
    }

    // We unclass, starting at 1
    k = 1;
    x_unclass[0] = k;

    double xi, xim1 = x_read[0];
    x_unik.push_back(xim1 + x_min);

    for (int i = 1; i < n; ++i)
    {
        xi = x_read[i];
        if (xi != xim1)
        {
            ++k;
            x_unik.push_back(static_cast<double>(xi + x_min));
        }
        x_unclass[i] = k;
        xim1 = xi;
    }

    // We put into the right order
    for (int i = 0; i < n; ++i)
    {
        x_uf[x_order[i]] = x_unclass[i];
    }
}

void quf_int(int n, int *x_uf, void *px, vector<double> &x_unik, int x_min, int max_value, bool is_double = false)
{
    // Principle:
    // we go through x only once
    // we keep a table of all x values
    // whenever we find a new value of x, we save it
    // px: pointer to the values of x -- R vector READ ONLY!!!
    // px: either a R vector of type integer // either a R vector of type double
    // px: we know it thanks to the value of 'is_double'

    int *px_int = (int *)px;
    double *px_dble = (double *)px;

    int n_unik = 0; // nber of uniques minus one

    // radix array
    vector<int> x_lookup(max_value + 1, 0);

    int x_tmp, x_pos;
    for (int i = 0; i < n; ++i)
    {
        x_tmp = is_double ? static_cast<int>(px_dble[i]) - x_min : px_int[i] - x_min;

        if (x_lookup[x_tmp] == 0)
        {
            ++n_unik;
            x_uf[i] = n_unik;
            x_unik.push_back(is_double ? px_dble[i] : static_cast<double>(px_int[i]));
            x_lookup[x_tmp] = n_unik;
        }
        else
        {
            x_pos = x_lookup[x_tmp];
            x_uf[i] = x_pos;
        }
    }
}

[[cpp11::register]] list cpp_quf_gnl(SEXP x)
{

    // INT: we try as possible to send the data to quf_int, the most efficient function
    // for data of large range, we have a separate algorithms that avoids the creation
    // of a possibly too large lookup table.
    // for extreme ranges (> 2**31) we use the algorithm for double

    // DOUBLE: quf_int is the most efficient function.
    // even if we have a vector of double, we check whether the underlying structure
    // is in fact int, so that we can send it to quf_int
    // if the vector is really double, then the checking cost is negligible
    // if the vector is in fact int, the cost of checking is largely compensated by the
    //   efficiency of the algorithm

    // STRING: for string vectors, we first transform them into ULL using their pointer
    // before applying them the algorithm for doubles
    // note that the x_unik returned is NOT of type string but is equal to
    // the location of the unique elements

    int n = Rf_length(x);

    // vector<int> x_uf(n);
    SEXP r_x_uf = PROTECT(Rf_allocVector(INTSXP, n));
    int *x_uf = INTEGER(r_x_uf);
    vector<double> x_unik;

    // preparation for strings
    bool IS_STR = false;
    vector<unsigned long long> x_ull;

    bool IS_INT = false;
    bool is_int_in_double = false;
    if (TYPEOF(x) == REALSXP)
    {
        // we check if underlying structure is int
        IS_INT = true;
        double *px = REAL(x);
        for (int i = 0; i < n; ++i)
        {
            if (!(px[i] == (int)px[i]))
            {
                IS_INT = false;
                break;
            }
        }

        is_int_in_double = IS_INT; // true: only if x is REAL + INT test OK
    }
    else if (TYPEOF(x) == STRSXP)
    {

        IS_STR = true;
        uintptr_t xi_uintptr;

        for (int i = 0; i < n; ++i)
        {
            const char *pxi = CHAR(STRING_ELT(x, i));
            xi_uintptr = reinterpret_cast<uintptr_t>(pxi);

            x_ull.push_back(static_cast<unsigned long long>(xi_uintptr));

            // Rcout << xi_uintptr << "  ----  " << xi_ull << "\n";
        }
    }
    else
    {
        IS_INT = true;
    }

    if (IS_INT)
    {
        // integer

        double max_value;
        void *px_generic;
        int X_MIN;
        if (is_int_in_double)
        {
            double *px = REAL(x);
            px_generic = REAL(x);
            double x_min = px[0], x_max = px[0], x_tmp;
            for (int i = 1; i < n; ++i)
            {
                x_tmp = px[i];
                if (x_tmp > x_max)
                    x_max = x_tmp;
                if (x_tmp < x_min)
                    x_min = x_tmp;
            }
            X_MIN = static_cast<int>(x_min);
            max_value = x_max - x_min;
        }
        else
        {
            int *px = INTEGER(x);
            px_generic = INTEGER(x);
            int x_min = px[0], x_max = px[0], x_tmp;
            for (int i = 1; i < n; ++i)
            {
                x_tmp = px[i];
                if (x_tmp > x_max)
                    x_max = x_tmp;
                if (x_tmp < x_min)
                    x_min = x_tmp;
            }
            X_MIN = x_min;
            max_value = x_max - x_min;
        }

        // creating + copying a 10**5 vector takes about 0.5ms which is OK even if n << 10**5
        // quf_int is really extremely efficient, so we use it whenever appropriate
        // In quf_int_gnl we create two n-size vectors + the method is less efficient by construction
        // so we go into quf_int whenever max_value <= 2.5*n

        if (max_value < 100000 || max_value <= 2.5 * n)
        {
            quf_int(n, x_uf, px_generic, x_unik, X_MIN, static_cast<int>(max_value), is_int_in_double);
        }
        else if (max_value < 0x10000000)
        {
            // we don't cover ranges > 2**31 (uints are pain in the neck)
            quf_int_gnl(n, x_uf, px_generic, x_unik, X_MIN, is_int_in_double);
        }
        else
        {
            // ranges > 2**31 => as double

            if (is_int_in_double)
            {
                quf_double(n, x_uf, (double *)px_generic, x_unik);
            }
            else
            {
                // we need to create a vector of double, otherwise: pointer issue
                vector<double> x_dble(n);
                int *px = INTEGER(x);
                for (int i = 0; i < n; ++i)
                    x_dble[i] = static_cast<double>(px[i]);
                quf_double(n, x_uf, x_dble.data(), x_unik);
            }
        }
    }
    else if (IS_STR)
    {
        // string -- beforehand transformed as ULL
        quf_double(n, x_uf, x_ull.data(), x_unik, true);
    }
    else
    {
        // double
        double *px = REAL(x);
        quf_double(n, x_uf, px, x_unik);
    }

    UNPROTECT(1);

    return list({"x_uf"_nm = r_x_uf,
                 "x_unik"_nm = x_unik});
}

void quf_single(void *px_in, string &x_type, int n, int *x_uf, vector<double> &x_unik)
{

    // preparation for strings
    bool IS_STR = x_type == "string";

    bool IS_INT = false;
    bool is_int_in_double = false;
    if (x_type == "double")
    {
        // we check if underlying structure is int
        IS_INT = true;
        double *px = (double *)px_in;
        for (int i = 0; i < n; ++i)
        {
            if (!(px[i] == (int)px[i]))
            {
                IS_INT = false;
                break;
            }
        }

        is_int_in_double = IS_INT; // true: only if x is REAL + INT test OK
    }
    else if (x_type == "int")
    {
        IS_INT = true;
    }

    if (IS_INT)
    {
        // integer

        double max_value;
        void *px_generic;
        int X_MIN;
        if (is_int_in_double)
        {
            double *px = (double *)px_in;
            px_generic = (double *)px_in;
            double x_min = px[0], x_max = px[0], x_tmp;
            for (int i = 1; i < n; ++i)
            {
                x_tmp = px[i];
                if (x_tmp > x_max)
                    x_max = x_tmp;
                if (x_tmp < x_min)
                    x_min = x_tmp;
            }
            X_MIN = static_cast<int>(x_min);
            max_value = x_max - x_min;
        }
        else
        {
            int *px = (int *)px_in;
            px_generic = (int *)px_in;
            int x_min = px[0], x_max = px[0], x_tmp;
            for (int i = 1; i < n; ++i)
            {
                x_tmp = px[i];
                if (x_tmp > x_max)
                    x_max = x_tmp;
                if (x_tmp < x_min)
                    x_min = x_tmp;
            }
            X_MIN = x_min;
            max_value = x_max - x_min;
        }

        // creating + copying a 10**5 vector takes about 0.5ms which is OK even if n << 10**5
        // quf_int is really extremely efficient, so we use it whenever appropriate
        // In quf_int_gnl we create two n-size vectors + the method is less efficient by construction
        // so we go into quf_int whenever max_value <= 2.5*n

        if (max_value < 100000 || max_value <= 2.5 * n)
        {
            quf_int(n, x_uf, px_generic, x_unik, X_MIN, static_cast<int>(max_value), is_int_in_double);
        }
        else if (max_value < 0x10000000)
        {
            // we don't cover ranges > 2**31 (uints are pain in the neck)
            quf_int_gnl(n, x_uf, px_generic, x_unik, X_MIN, is_int_in_double);
        }
        else
        {
            // ranges > 2**31 => as double

            if (is_int_in_double)
            {
                quf_double(n, x_uf, (double *)px_generic, x_unik);
            }
            else
            {
                // we need to create a vector of double, otherwise: pointer issue
                vector<double> x_dble(n);
                int *px = (int *)px_in;
                for (int i = 0; i < n; ++i)
                    x_dble[i] = static_cast<double>(px[i]);
                quf_double(n, x_uf, x_dble.data(), x_unik);
            }
        }
    }
    else if (IS_STR)
    {
        // string -- beforehand transformed as ULL
        quf_double(n, x_uf, px_in, x_unik, true);
    }
    else
    {
        // double
        quf_double(n, x_uf, px_in, x_unik);
    }
}

