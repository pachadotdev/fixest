#include "08_0_strings.hpp"

// we change a:b into a*b
// we don't touch what's in paren: I(a:b) stays I(a:b)
// x = c("x1:x2:a(6:7)", "x5", "i(aa, 5:6):jjl", "base::poly(x, 5)")
string colon_to_star_single(const char *str)
{
    string res = "";

    int n = strlen(str);

    int i = 0;
    bool in_quote = false;
    char quote = '"';
    int n_paren = 0;
    while (i < n)
    {

        if (in_quote)
        {
            if (str[i] == quote)
            {
                in_quote = false;
            }
        }
        else if (str[i] == '"' || str[i] == '\'')
        {
            in_quote = true;
            quote = str[i];
        }
        else if (n_paren > 0)
        {
            if (str[i] == '(')
            {
                ++n_paren;
            }
            else if (str[i] == ')')
            {
                --n_paren;
            }
        }
        else if (str[i] == '(')
        {
            ++n_paren;
        }
        else if (str[i] == ':')
        {
            if (i + 1 < n && str[i + 1] != ':' && i - 1 >= 0 && str[i - 1] != ':')
            {
                // OK
                res += '*';
                ++i;
                continue;
            }
        }

        if (i == n)
            stop("Error in the index.");

        res += str[i++];
    }

    return res;
}

[[cpp11::register]] strings cpp_colon_to_star_(SEXP Rstr)
{

    int n = LENGTH(Rstr);

    writable::strings res(n);
    for (int i = 0; i < n; ++i)
    {
        const char *str = CHAR(STRING_ELT(Rstr, i));
        res[i] = colon_to_star_single(str);
    }

    return res;
}
