# small utilities ----

escape_regex = function(x){
    # escape special characters in regular expressions => to make it as "fixed"

    res = gsub("((?<=[^\\\\])|(?<=^))(\\$|\\.|\\+|\\(|\\)|\\[|\\]|\\?|\\^)", "\\\\\\2", x, perl = TRUE)
    res
}

.addCommas = function(x){

    if(!is.finite(x)) return(as.character(x))

    cpp_add_commas(x, r = 1, whole = TRUE)
}

addCommas = function(x){
    sapply(x, .addCommas)
}

decimalFormat = function(x){

    who_valid = which(!is.na(x) & is.numeric(x))
    if(length(who_valid) == 0) return(x)

    res = x

    x_valid = x[who_valid]
    xPower = log10(abs(x_valid))

    if(min(xPower) > 0){
        pow_round = max(1, 6 - ceiling(xPower))
    } else if(min(xPower) < -5){
        pow_round = ceiling(abs(min(xPower))) + 2
    } else {
        pow_round = 6
    }

    res[who_valid] = round(x_valid, pow_round)

    res
}

numberFormat_single = function(x, type = "normal"){
    # For numbers higher than 1e9 => we apply a specific formatting
    # idem for numbers lower than 1e-4

    if(is.character(x)) return(x)

    if(is.na(x)){
        return(NA)
    }

    if(x == 0) return("0")

    exponent = floor(log10(abs(x)))

    if(-4 < exponent && exponent < 9){
        if(exponent > 0){
            return(addCommas(x))
        } else {
            return(as.character(decimalFormat(x)))
        }

    }

    left_value = round(x*10**-exponent, 3)

    if(type == "latex"){
        res = paste0("$", left_value, "\\times 10^{", exponent, "}$")
    } else {
        res = paste0(left_value, "e", ifelse(exponent > 0, "+", ""), exponent)
    }

    res
}

numberFormatLatex = function(x){
    sapply(x, numberFormat_single, type = "latex")
}

numberFormatNormal = function(x){
    sapply(x, numberFormat_single)
}

mysignif = function (x, d = 2, r = 1){

    .mysignif = function(x, d, r) {
        if (is.na(x)) {
            return(NA)
        }

        if(abs(x) >= 10^(d - 1)){
            return(round(x, r))
        } else {
            return(signif(x, d))
        }
    }
    sapply(x, .mysignif, d = d, r = r)
}

format_nber_single = function(x, digits, round = FALSE, pow_above = 10, pow_below = -5, tex = FALSE){
    # Attention aux nombres ronds => pas chiffre apres la virgule!

    if(is.na(x) || !is.numeric(x)){
        return(x)
    }

    if(round){
        # super simple
        return(cpp_add_commas(x, digits))

    } else {
        whole = (x %% 1) == 0
        if(abs(x) >= 10^(digits - 1)){
            x = round(x, 1)
        } else {
            x = signif(x, digits)
        }
    }

    exponent = floor(log10(abs(x)))

    if(exponent >= pow_above || exponent <= pow_below){
        left_value = round(x*10**-exponent, min(digits, 2))

        if(tex){
            res = paste0("$", left_value, "\\times 10^{", exponent, "}$")
        } else {
            res = paste0(left_value, "e", ifelse(exponent > 0, "+", ""), exponent)
        }

    } else if(exponent >= 0){
        r = max(-(exponent + 1 - digits), 1)
        res = cpp_add_commas(x, r, whole)

    } else if(abs(x) > 10**(-digits)){
        res = sprintf("%.*f", digits, x)

    } else {
        res = sprintf("%.*f", abs(exponent), x)
    }

    res
}

format_number = function(x, digits = 4, round = FALSE, pow_above = 10, pow_below = -5, tex = FALSE){
    sapply(x, format_nber_single, digits = digits, round = round, pow_above = pow_above, pow_below = pow_below, tex = tex)
}

index_2D_to_1D = function(i, j, n_j) 1 + (i - 1) * n_j + j - 1

par_fit = function(my_par, id){
    # simple function that extends the plot parameters
    my_par = rep(my_par, ceiling(max(id) / length(my_par)))
    my_par[id]
}

dict_apply = function(x, dict = NULL){

    check_arg(dict, "NULL named character vector no na", .message = "The argument 'dict' must be a dictionnary, ie a named vector (eg dict=c(\"old_name\"=\"new_name\")")

    if(missing(dict) || length(dict) == 0){
        return(x)
    }

    # We make the dictionary names space agnostic, adds a handful of us only
    if(any(grepl(" ", x, fixed = TRUE))){
        x_tmp = gsub(" ", "", x, fixed = TRUE)
    } else {
        x_tmp = x
    }

    if(any(grepl(" ", names(dict), fixed = TRUE))){
        names(dict) = gsub(" ", "", names(dict), fixed = TRUE)
        if(anyDuplicated(names(dict))){
            dup = duplicated(names(dict))
            stop("The dictionary 'dict' contains several items with the same names, it concerns ",
                 enumerate_items(names(dict)[dup]), " (note that spaces are ignored).")
        }
    }

    who = x_tmp %in% names(dict)
    x[who] = dict[as.character(x_tmp[who])]
    x
}

keep_apply = function(x, keep = NULL, logical = FALSE){

    if(missing(keep) || length(keep) == 0){
        if(logical){
            return(rep(TRUE, length(x)))
        } else {
            return(x)
        }
    }

    check_arg(keep, "character vector no na", .message = "The arg. 'keep' must be a vector of regular expressions (see help(regex)).")

    res = x

    qui_keep = rep(FALSE, length(x))
    for(var2keep in keep){
        if(grepl("^!%", var2keep)) var2keep = gsub("^!%", "%!", var2keep)

        vect2check = res
        if(grepl("^%", var2keep)){
            var2keep = gsub("^%", "", var2keep)
            vect2check = names(res)
            if(is.null(vect2check)){
                # warning("In keep, the special character '%' cannot be used here.")
                vect2check = res
            }
        }

        if(grepl("^!", var2keep)){
            qui_keep = qui_keep | !grepl(substr(var2keep, 2, nchar(var2keep)), vect2check)
        } else {
            qui_keep = qui_keep | grepl(var2keep, vect2check)
        }
    }

    if(logical) return(qui_keep)

    res[qui_keep]
}

drop_apply = function(x, drop = NULL){

    if(missing(drop) || length(drop) == 0){
        return(x)
    }

    check_arg(drop, "character vector no na", .message = "The arg. 'drop' must be a vector of regular expressions (see help(regex)). ")

    res = x

    for(var2drop in drop){
        if(grepl("^!%", var2drop)) var2drop = gsub("^!%", "%!", var2drop)

        vect2check = res
        if(grepl("^%", var2drop)){
            var2drop = gsub("^%", "", var2drop)
            vect2check = names(res)
            if(is.null(vect2check)){
                # warning("In drop, the special character '%' cannot be used here.")
                vect2check = res
            }
        }

        if(grepl("^!", var2drop)){
            res = res[grepl(substr(var2drop, 2, nchar(var2drop)), vect2check)]
        } else {
            res = res[!grepl(var2drop, vect2check)]
        }
    }

    res
}

order_apply = function(x, order = NULL){

    if(missing(order) || length(order) == 0){
        return(x)
    }

    check_arg(order, "character vector no na", .message = "The arg. 'order' must be a vector of regular expressions (see help(regex)). ")

    res = x

    for(var2order in rev(order)){
        if(grepl("^!%", var2order)) var2order = gsub("^!%", "%!", var2order)

        vect2check = res
        if(grepl("^%", var2order)){
            var2order = gsub("^%", "", var2order)
            vect2check = names(res)
            if(is.null(vect2check)){
                # warning("In order, the special character '%' cannot be used here.")
                vect2check = res
            }
        }

        if(grepl("^!", var2order)){
            who = !grepl(substr(var2order, 2, nchar(var2order)), vect2check)
            res = c(res[who], res[!who])
        } else {
            who = grepl(var2order, vect2check)
            res = c(res[who], res[!who])
        }
    }

    res
}

charShorten = function(x, width){
    # transforms characters => they can't go beyond a certain width
    # two dots are added to suggest longer character
    # charShorten("bonjour", 5) => "bon.."
    n = nchar(x)

    if(n > width && width > 3){
        res = substr(x, 1, width - 2)
        res = paste0(res, "..")
    } else {
        res = x
    }

    res
}


show_vars_limited_width = function(charVect, nbChars = 60, addS = FALSE){
    # There are different cases

    n = length(charVect)

    if(n == 1){
        return(charVect)
    }

    text_s = ifelse(addS, "s ", "")

    nb_char_vect = nchar(charVect)
    sumChar = cumsum(nb_char_vect) + (0:(n-1))*2 + 3 + 1

    if(max(sumChar) < nbChars){
        text = paste0(text_s, paste0(charVect[-n], collapse = ", "), " and ", charVect[n])
        return(text)
    }

    qui = max(which.max(sumChar > nbChars - 8) - 1, 1)

    nb_left = n - qui

    text = paste0(text_s, paste0(charVect[1:qui], collapse = ", "), " and ", nb_left, " other", ifelse(nb_left>1, "s", ""), ".")

    return(text)
}


char2num = function(x, addItem = FALSE){
    # we transform the data to numeric => faster analysis

    # special case
    qui = which(x == "")
    if(length(qui) > 0){
        x[qui] = "xxEMPTYxx"
    }

    x_unik = unique(x)
    dict = 1:length(x_unik)
    names(dict) = x_unik
    x_num = dict[x]

    names(x_num) = NULL

    if(addItem){
        res = list(x = x_num, items = x_unik)
        return(res)
    } else {
        return(x_num)
    }

}

quickUnclassFactor = function(x, addItem = FALSE, sorted = FALSE){
    # does as unclass(as.factor(x))
    # but waaaaay quicker

    not_num = !is.numeric(x)
    is_char_convert = not_num && !is.character(x)

    if(is_char_convert){
        res = cpp_quf_gnl(as.character(x))
    } else {
        res = cpp_quf_gnl(x)
    }

    if(sorted){

        if(not_num){
            items = x[res$x_unik]
        } else {
            items = res$x_unik
        }

        x = res$x_uf

        new_order = order(items)
        order_new_order = order(new_order)
        x_uf = order_new_order[x]

        if(addItem){
            if(is_char_convert){
                res = list(x = x_uf, items = as.character(items[new_order]))
            } else {
                res = list(x = x_uf, items = items[new_order])
            }

            return(res)
        } else {
            return(x_uf)
        }
    }

    if(addItem){

        if(not_num){
            if(is_char_convert){
                items = as.character(x[res$x_unik])
            } else {
                items = x[res$x_unik]
            }

            res = list(x = res$x_uf, items = items)
        } else {
            names(res) = c("x", "items")
        }

        return(res)
    } else {
        return(res$x_uf)
    }
}

missnull = function(x) missing(x) || is.null(x)

MISSNULL = function(x){
    # same a missnull but we also check evaluation

    if(missing(x)) return(TRUE)

    # we check evaluation and nullity
    error_sender(x, up = 1, arg_name = deparse(substitute(x)))

    is.null(x)
}

isScalar = function(x, int = FALSE) {
    if(length(x) == 1L && is.numeric(x) && is.finite(x)){
        if(int){
            return(x %% 1 == 0)
        } else {
            return(TRUE)
        }
    } else {
        return(FALSE)
    }
}

isLogical = function(x) length(x) == 1L && is.logical(x) && !is.na(x)

isSingleChar = function(x) length(x) == 1L && is.character(x) && !is.na(x)

fml2varnames = function(fml, combine_fun = FALSE){
    # This function transforms a one sided formula into a
    # character vector for each variable

    # In theory, I could just use terms.formula to extract the variable.
    # But I can't!!!! Because of this damn ^ argument.
    # I need to apply a trick

    # combine_fun: whether to add a call to combine_clusters_fast

    # Only the ^ users "pay the price"

    if("^" %in% all.vars(fml, functions = TRUE)){
        # new algo using fml_combine, more robust
        fml_char = as.character(fml)[2]
        all_var_names = fml_combine(fml_char, TRUE, vars = TRUE)
        if(!combine_fun){
            all_var_names = rename_hat(all_var_names)
        }

    } else {
        t = terms(fml)
        all_var_names = attr(t, "term.labels")
        # => maybe I should be more cautious below???? Sometimes ':' really means stuff like 1:5
        all_var_names = gsub(":", "*", all_var_names) # for very special cases
    }


    all_var_names
}

msg_na_inf = function(any_na, any_inf){
    if(any_na && any_inf){
        res = "NA and infinite values"
    } else if(any_na){
        res = "NA values"
    } else {
        res = "infinite values"
    }

    res
}

extract_fe_slope = function(t){
    # input: x$fixef_terms which contains the slopes and the fixef

    fixef_vars = t[!grepl("\\[", t)]
    slope_terms = t[grepl("\\[", t)]
    slope_vars = gsub(".+\\[|\\]", "", slope_terms)
    slope_fe = gsub("\\[.+", "", slope_terms)
    fe_all = gsub("\\[.+", "", t)

    list(fixef_vars=fixef_vars, slope_vars=slope_vars, slope_fe=slope_fe, slope_terms=slope_terms, fe_all=fe_all)
}

deparse_long = function(x){
    dep_x = deparse(x, width.cutoff = 500)
    if(length(dep_x) == 1){
        return(dep_x)
    } else {
        return(paste(gsub("^ +", "", dep_x), collapse = ""))
    }
}

isVector = function(x){
    # it seems that when you subselect in data.table
    # sometimes it does not yield a vector
    # so i cannot use is.vector to check the consistency

    if(is.vector(x)){
        return(TRUE)
    } else {
        if(is.null(dim(x)) && !is.list(x)){
            return(TRUE)
        }
    }
    return(FALSE)
}


fill_with_na = function(x, object){
    if(is.null(object$obs_selection)){
        return(x)
    }

    res = rep(NA, object$nobs_origin)
    qui = 1:object$nobs_origin

    for(i in seq_along(object$obs_selection)){
        qui = qui[object$obs_selection[[i]]]
    }

    res[qui] = x

    return(res)
}

trim_obs_removed = function(x, object){
    if(is.null(object$obs_selection)){
        return(x)
    }

    for(i in seq_along(object$obs_selection)){
        x = x[object$obs_selection[[i]]]
    }

    x
}

is_operator = function(x, op) if(length(x) <= 1) FALSE else x[[1]] == op

fml_breaker = function(fml, op){
    res = list()
    k = 1
    while(is_operator(fml, op)){
        res[[k]] = fml[[3]]
        k = k + 1
        fml = fml[[2]]
    }
    res[[k]] = fml

    res
}

fml_maker = function(lhs, rhs){

    while(is_operator(lhs, "(")){
        lhs = lhs[[2]]
    }

    if(missing(rhs)){
        if(is_operator(lhs, "~")){
            return(lhs)
        }
        res = ~ .
        res[[2]] = lhs
    } else {

        while(is_operator(rhs, "(")){
            rhs = rhs[[2]]
        }

        res = . ~ .
        res[[2]] = lhs
        res[[3]] = rhs
    }

    res
}

fml_split_internal = function(fml, split.lhs = FALSE){

    fml_split_tilde = fml_breaker(fml, "~")
    k = length(fml_split_tilde)

    # NOTA in fml_breaker: to avoid copies, the order of elements returned is reversed

    # currently res is the LHS
    res = list(fml_split_tilde[[k]])

    if(k == 2){
        rhs = fml_breaker(fml_split_tilde[[1]], "|")
        l = length(rhs)

    } else if(k == 3){
        rhs  = fml_breaker(fml_split_tilde[[2]], "|")
        l = length(rhs)
        rhs_right = fml_breaker(fml_split_tilde[[1]], "|")

        if(length(rhs_right) > 1){
            stop_up("Problem in the formula: the formula in the RHS (expressing the IVs) cannot be multipart.")
        }

        # The rightmost element of the RHS is in position 1!!!!
        iv_fml = fml_maker(rhs[[1]], rhs_right[[1]])

        rhs[[1]] = iv_fml

    } else {
        # This is an error
        stop_up("Problem in the formula: you cannot have more than one RHS part containing a formula.")
    }

    if(!split.lhs){
        new_fml = fml_maker(res[[1]], rhs[[l]])

        res[[1]] = new_fml
        if(l == 1) return(res)
        res[2:l] = rhs[(l - 1):1]

    } else {
        res[1 + (1:l)] = rhs[l:1]
    }

    res
}


fml_split = function(fml, i, split.lhs = FALSE, text = FALSE, raw = FALSE){
    # I had to create that fonction to cope with the following formula:
    #
    #                       y ~ x1 | fe1 | u ~ z
    #
    # For Formula to work well, one would need to write insstead: y | x1 | fe1 | (u ~ z)
    # that set of parentheses are superfluous

    my_split = fml_split_internal(fml, split.lhs)

    if(raw){
        return(my_split)
    } else if(text){

        if(!missing(i)){
            return(deparse_long(my_split[[i]]))
        } else {
            return(sapply(my_split, deparse_long))
        }

    } else if(!missing(i)) {
        return(fml_maker(my_split[[i]]))

    } else {
        res = lapply(my_split, fml_maker)

        return(res)
    }

}

error_sender = function(expr, ..., clean, up = 0, arg_name){
    res = tryCatch(expr, error = function(e) structure(list(conditionCall(e), conditionMessage(e)), class = "try-error"))

    if("try-error" %in% class(res)){
        set_up(1 + up)
        msg = paste0(..., collapse = "")

        if(nchar(msg) == 0){
            if(missing(arg_name)){
                arg_name = deparse(substitute(expr))
            }
            msg = paste0("Argument '", arg_name, "' could not be evaluated: ")
            stop_up(msg, res[[2]])

        } else {

            # The expression that is evaluated is absolutely non informative for the user
            call_non_informative = deparse(substitute(expr), 100)[1]
            call_error = deparse(res[[1]], 100)[1]

            if(call_error == call_non_informative || call_error == "NULL" ||
               grepl("^(doTry|eval)", call_error)){
                call_error = ""

            } else {
                call_error = paste0("In ", call_error, ": ")
            }

            err = res[[2]]

            if(grepl("^in eval\\(str[^:]+:\n", err)){
                err = sub("^in eval\\(str[^:]+:\n", "", err)
            }

            if(!missing(clean)){

                if(grepl(" => ", clean)){
                    clean_split = strsplit(clean, " => ")[[1]]
                    from = clean_split[1]
                    to = clean_split[2]
                } else {
                    from = clean
                    to = ""
                }

                stop_up(msg, "\n  ", call_error, gsub(from, to, err))
            } else {
                stop_up(msg, "\n  ", call_error, err)
            }
        }
    }

    res
}

is_fml_inside = function(fml){
    # we remove parentheses first

    while(is_operator(fml, "(")){
        fml = fml[[2]]
    }

    is_operator(fml, "~")
}


merge_fml = function(fml_linear, fml_fixef = NULL, fml_iv = NULL){

    is_fe = length(fml_fixef) > 0
    is_iv = length(fml_iv) > 0

    if(!is_fe && !is_iv){
        res = fml_linear
    } else {
        fml_all = deparse_long(fml_linear)

        if(is_fe){
            # we add parentheses if necessary
            if(is_operator(fml_fixef[[2]], "|")){
                fml_all[[2]] = paste0("(", as.character(fml_fixef)[2], ")")
            } else {
                fml_all[[2]] = as.character(fml_fixef)[2]
            }
        }

        if(is_iv) fml_all[[length(fml_all) + 1]] = deparse_long(fml_iv)

        res = as.formula(paste(fml_all, collapse = "|"))
    }

    res
}


fixest_fml_rewriter = function(fml){
    # Currently performs the following
    # - expands lags
    # - protects powers: x^3 => I(x^3)
    #
    # fml = sw(f(y, 1:2)) ~ x1 + l(x2, 1:2) + x2^2 | fe1 | y ~ z::e + g^3
    # fml = y ~ 1 | id + period | l(x_endo, -1:1) ~ l(x_exo, -1:1)

    fml_text = deparse_long(fml)

    isPanel = grepl("(^|[^\\._[:alnum:]])(f|d|l)\\(", fml_text)
    isPower = grepl("^", fml_text, fixed = TRUE)

    if(isPanel){
        # We rewrite term-wise

        fml_parts = fml_split(fml, raw = TRUE)
        n_parts = length(fml_parts)

        #
        # LHS
        #

        # only panel: no power (bc no need), no interact

        # We tolerate multiple LHS and expansion
        lhs_text = fml_split(fml, 1, text = TRUE, split.lhs = TRUE)

        if(grepl("^(c|c?sw0?|list)\\(", lhs_text)){
            lhs_text2eval = gsub("^(c|c?sw0?|list)\\(", "sw(", lhs_text)
            lhs_names = eval(str2lang(lhs_text2eval))
        } else {
            lhs_names = lhs_text
        }

        lhs_all = error_sender(expand_lags_internal(lhs_names),
                               "Problem in the formula regarding lag/leads: ", clean = "__expand")

        if(length(lhs_all) > 1){
            lhs_fml = paste("c(", paste(lhs_all, collapse = ","), ")")
        } else {
            lhs_fml = lhs_all
        }

        lhs_text = lhs_fml

        #
        # RHS
        #

        # power + panel + interact

        if(isPower){
            # rhs actually also contains the LHS
            rhs_text = deparse_long(fml_parts[[1]])
            rhs_text = gsub("(?<!I\\()(\\b(\\.[[:alpha:]]|[[:alpha:]])[[:alnum:]\\._]*\\^[[:digit:]]+)", "I(\\1)", rhs_text, perl = TRUE)

            if(grepl("\\^[[:alpha:]]", rhs_text)){
                stop_up("The operator '^' between variables can be used only in the fixed-effects part of the formula. Otherwise, please use ':' instead.")
            }

            fml_rhs = as.formula(rhs_text)
        } else {
            fml_rhs = fml_maker(fml_parts[[1]])
        }

        rhs_terms = get_vars(fml_rhs)

        if(length(rhs_terms) == 0){
            rhs_text = "1"
        } else {
            rhs_terms = error_sender(expand_lags_internal(rhs_terms),
                                     "Problem in the formula regarding lag/leads: ", clean = "__expand")

            if(attr(terms(fml_rhs), "intercept") == 0){
                rhs_terms = c("-1", rhs_terms)
            }

            rhs_text = paste(rhs_terms, collapse = "+")
        }

        fml_linear = as.formula(paste0(lhs_text, "~", rhs_text))

        #
        # FE + IV
        #

        fml_fixef = fml_iv = NULL

        if(n_parts > 1){

            #
            # FE
            #

            # Only isPanel (although odd....)

            is_fe = !is_fml_inside(fml_parts[[2]])
            if(is_fe){

                if(identical(fml_parts[[2]], 0)){
                    fml_fixef = NULL

                } else {

                    fml_fixef = fml_maker(fml_parts[[2]])
                    fml_fixef_text = deparse_long(fml_fixef)

                    if(grepl("(l|d|f)\\(", fml_fixef_text)){
                        # We need to make changes
                        # 1st: for terms to work => we change ^ if present (sigh)

                        do_sub = grepl("^", fml_fixef_text, fixed = TRUE)

                        if(do_sub){
                            fml_fixef = as.formula(gsub("^", "%^%", fml_fixef_text, fixed = TRUE))
                        }

                        fixef_terms = attr(terms(fml_fixef), "term.labels")
                        fixef_text = error_sender(expand_lags_internal(fixef_terms),
                                                  "Problem in the formula regarding lag/leads: ", clean = "__expand")

                        if(do_sub){
                            fixef_text = gsub(" %^% ", "^", fixef_text, fixed = TRUE)
                        }

                        fml_fixef = as.formula(paste("~", paste(fixef_text, collapse = "+")))

                    }
                }
            }

            #
            # IV
            #

            if(n_parts == 3 || !is_fe){
                fml_iv = fml_maker(fml_parts[[n_parts]])

                fml_endo = .xpd(lhs = ~y, rhs = fml_iv[[2]])
                fml_inst = .xpd(lhs = ~y, rhs = fml_iv[[3]])

                endo_lag_expand = fixest_fml_rewriter(fml_endo)$fml
                inst_lag_expand = fixest_fml_rewriter(fml_inst)$fml

                fml_iv = .xpd(lhs = endo_lag_expand[[3]], rhs = inst_lag_expand[[3]])
            }
        }

        fml_new = merge_fml(fml_linear, fml_fixef, fml_iv)

    } else if(isPower){
        # This is damn costly.... but I have to deal with corner cases....
        # I think I should do that at the C level, this will be much faster I guess

        # We only take care of the RHS (we don't care about the LHS)
        no_lhs_text = gsub("^[^~]+~", "", fml_text)
        no_lhs_text = gsub("(?<!I\\()(\\b(\\.[[:alpha:]]|[[:alpha:]])[[:alnum:]\\._]*\\^[[:digit:]]+)", "I(\\1)", no_lhs_text, perl = TRUE)

        if(grepl("\\^[[:alpha:]]", no_lhs_text)){
            # We check if there is one ^ specifically in the RHS or in the IV part
            fml_parts = fml_split(fml, raw = TRUE)
            n_parts = length(fml_parts)

            rhs_txt = fml_split(fml, i = 1, text = TRUE)

            if(grepl("\\^[[:alpha:]]", rhs_txt)){
                stop_up("The operator '^' between variables can be used only in the fixed-effects part of the formula. Otherwise, please use ':' instead.")
            }

            if(is_fml_inside(fml_parts[[n_parts]])){
                iv_txt = fml_split(fml, i = n_parts, text = TRUE)

                if(grepl("\\^[[:alpha:]]", iv_txt)){
                    stop_up("The operator '^' between variables can be used only in the fixed-effects part of the formula. Otherwise, please use ':' instead.")
                }
            }

        }

        fml_new = as.formula(paste0(gsub("~.+", "", fml_text), "~", no_lhs_text))

    } else {
        res = list(fml = fml, isPanel = FALSE)
        return(res)
    }

    res = list(fml = fml_new, isPanel = isPanel)

    return(res)
}


check_set_types = function(x, types, msg){
    arg_name = deparse(substitute(x))
    check_arg(x, "os formula | character vector no na", .arg_name = arg_name, .up = 1)

    if("formula" %in% class(x)){
        x = attr(terms(x), "term.labels")
    }

    check_value_plus(x, "multi match", .choices = types, .arg_name = arg_name, .up = 1)

    x
}

check_set_digits = function(digits, up = 1){
    # Note that the argument name can be either digits or digits.stats

    set_up(up)
    check_value(digits, "integer scalar GE{1} | character scalar", .arg_name = deparse(substitute(digits)))

    if(is.character(digits)){
        d_type = substr(digits, 1, 1)
        d_value = substr(digits, 2, nchar(digits))

        # Control
        if(!d_type %in% c("r", "s")){
            arg = deparse(substitute(digits))
            stop_up("The argument '", arg, "' must start with 'r' (for round) or 's' (for significant). Currently it starts with '", d_type,"' which is not valid.\nExample of valid use: digits = 'r3'.")
        }

        round = d_type == "r"

        if(!grepl("^[0-9]$", d_value)){
            arg = deparse(substitute(digits))
            stop_up("The argument '", arg, "' must be equal to the character 'r' or 's' followed by a single digit. Currently '", digits,"' is not valid.\nExample of valid use: digits = '", d_type, "3'.")
        }

        digits = as.numeric(d_value)

    } else {
        round = FALSE
    }

    list(digits = digits, round = round)
}

get_vars = function(x){
    attr(terms(x), "term.labels")
}

mat_posdef_fix = function(X, tol = 1e-10){
    # X must be a symmetric matrix
    # We don't check it

    if(any(diag(X) < tol)){
        e = eigen(X)
        dm = dimnames(X)
        X = tcrossprod(e$vectors %*% diag(pmax(e$values, tol), nrow(X)), e$vectors)
        dimnames(X) = dm
    }

    return(X)
}


is_fixest_call = function(){
    nf = sys.nframe()

    if(nf < 5) return(FALSE)

    last_calls = sapply(tail(sys.calls(), 13), function(x) deparse(x)[1])

    any(grepl("fixest", last_calls[-length(last_calls)]))
}

all_vars_with_i_prefix = function(fml){
    # fml = a ~ x1^x2 + i(x3, i.x4) + x5*i(x6, I(i.x7))

    vars = all.vars(fml)
    if(any(grepl("^i\\..+", vars))){
        fml_dp = deparse_long(fml)
        # for i. to work it MUST be the second argument (var)
        # valid cases:
        # - i(x1, i.x2)
        # - i(var = i.x2, x1)
        # - i(x1, i.I(x7))
        # Not valid:
        # - i(x1, I(i.x7))
        #
        # This means that in the parsed formula, i. is always preceded by a space and
        # either a "," or a "="

        # Maybe later: add nice error messages reminding how to use i.

        qui_i = which(grepl("^i\\..+", vars))
        i_vars = vars[qui_i]
        for(i in seq_along(i_vars)){
            fml_split = strsplit(fml_dp, i_vars[i], fixed = TRUE)[[1]]
            n = length(fml_split) - 1
            for(j in 1:n){
                part = fml_split[j]
                if(grepl("(,|var =) *$", part)){
                    part = gsub("\\([^\\)]+\\)", "", part)
                    if(grepl("i\\(", part)){
                        # OK!
                        ii = qui_i[i]
                        vars[ii] = substr(vars[ii], 3, nchar(vars[ii]))
                    }
                }
            }
        }
    }

    vars
}


# function to normalize character vectors into variable names
as_varname = function(x){
    # "x1" => "x1"
    # "(x1)" => "`(x1)`"


    qui_pblm = grepl("[^[:alnum:]\\._]", x)
    if(any(qui_pblm)){
        x[qui_pblm] = paste0("`", x[qui_pblm], "`")
        x[qui_pblm] = gsub("``", "`", x[qui_pblm])
    }

    x
}

update_file = function(path, text){

    if(file.exists(path)){
        text_1 = readLines(path)

        text_clean_1 = unlist(strsplit(text_1, "\n"))
        text_clean_2 = unlist(strsplit(text, "\n"))

        text_clean_1 = text_clean_1[grepl("[[:alnum:][:punct:]]", text_clean_1)]
        text_clean_2 = text_clean_2[grepl("[[:alnum:][:punct:]]", text_clean_2)]

        do_write = length(text_clean_1) != length(text_clean_2) || any(text_clean_1 != text_clean_2)

    } else {
        do_write = TRUE
    }

    # We write only if the text is different
    if(do_write){
        message("writing '", path, "'", appendLF = FALSE)
        f = file(path, "w", encoding = "utf-8")
        writeLines(text, f)
        close(f)
        message(".")
    }
}


fetch_arg_deparse = function(arg){
    # Utility to find out what the user originally provided as argument
    # Only useful for functions that are deeply nested (typically vcov)

    sc = rev(sys.calls())

    # initialization of the name
    arg_name = deparse_long(sc[[2]][[arg]])

    n = length(sc)
    if(n > 2){
        for(i in 2:length(sc)){
            mc = sc[[i]]
            if(grepl("[tT]ry", mc[[1]])) next
            if(!arg %in% names(mc)) break
            arg_name = deparse_long(mc[[arg]])
        }
    }

    arg_name
}



# Simple utilities to avoid writing too much
# adds a small overhead but that's OK (about 1us which is OK)
any_missing = function(x1, x2, x3, x4, x5, x6){
    n = length(sys.call()) - 1
    switch(n,
           "1" = missing(x1),
           "2" = missing(x1) || missing(x2),
           "3" = missing(x1) || missing(x2) || missing(x3),
           "4" = missing(x1) || missing(x2) || missing(x3) || missing(x4),
           "5" = missing(x1) || missing(x2) || missing(x3) || missing(x4) || missing(x5),
           "6" = missing(x1) || missing(x2) || missing(x3) || missing(x4) || missing(x5) || missing(x6))
}

all_missing = function(x1, x2, x3, x4, x5, x6){
    n = length(sys.call()) - 1
    switch(n,
           "1" = missing(x1),
           "2" = missing(x1) && missing(x2),
           "3" = missing(x1) && missing(x2) && missing(x3),
           "4" = missing(x1) && missing(x2) && missing(x3) && missing(x4),
           "5" = missing(x1) && missing(x2) && missing(x3) && missing(x4) && missing(x5),
           "6" = missing(x1) && missing(x2) && missing(x3) && missing(x4) && missing(x5) && missing(x6))
}

use_t_distr = function(x){
    # whether to use the t-distribution or the normal
    # x: fixest estimation
    x$method %in% "feols" || (x$method %in% "feglm" && !x$family$family %in% c("poisson", "binomial"))
}


is_fun_in_char = function(fml_char, fun){
    extract_fun(fml_char, fun, bool = TRUE)
}

extract_fun = function(fml_char, fun, err_msg = NULL, bool = FALSE, drop = TRUE){
    # A formula, in character form // can be a vector
    # fun: a function name in regex form
    #
    # returns a list:
    # before
    # fun
    # after
    #
    # fml_char = "to_integer(id, sw0(fe))" ;  fun = "sw0?"

    only_one = TRUE

    regex = paste0(c("(?<=^)", "(?<=[^[:alnum:]\\._])"),
                   fun, "\\(", collapse = "|")

    is_there = grepl(regex, fml_char, perl = TRUE)

    if(bool) return(is_there)

    res = list()

    n = length(fml_char)
    for(i in 1:n){
        # fml_i = fml_char[1]

        fml_i = fml_char[i]
        if(is_there[i]){

            fml_split = strsplit(fml_i, regex, perl = TRUE)[[1]]

            if(only_one && length(fml_split) > 2){
                if(is.null(err_msg)){
                    stop_up("Only one function '", fun, "' can be used at a time.")
                } else {
                    stop_up(err_msg)
                }
            }

            # we need to add a last character otherwise => problems
            fml_value = paste0(fml_split[2], " ")
            fml_right_split = strsplit(fml_value, ")", fixed = TRUE)[[1]]
            n_open = pmax(lengths(strsplit(fml_right_split, "(", fixed = TRUE)) - 1, 0)

            n_parts = length(fml_right_split)
            fun_closed = which((1 + cumsum(n_open) - 1:n_parts) == 0)

            fun_name = substr(fml_i, nchar(fml_split[1]) + 1, nchar(fml_i))
            fun_name = gsub("\\(.+", "", fun_name)

            fun_char = paste0(fun_name, "(", paste(fml_right_split[1:fun_closed], collapse = ")"), ")")

            fml_right_rest = trimws(paste0(fml_right_split[-(1:fun_closed)], collapse = ")"))

            res[[i]] = list(before = fml_split[1],
                            fun = fun_char,
                            after = fml_right_rest)


        } else {
            res[[i]] = list(before = fml_i,
                            fun = "",
                            after = "")
        }
    }

    if(n == 1 && drop) res = res[[1]]

    return(res)
}

is_numeric_in_char = function(x){
    res = tryCatch(as.numeric(x), warning = function(x) "not numeric")
    !identical(res, "not numeric")
}

insert_in_between = function(x, y){
    n_x = length(x)
    n_y = length(y)

    if(n_y == 1) y = rep(y, n_x)
    if(n_x == 1) x = rep(x, n_y)
    n_x = length(x)

    res = rep(x, each = 2)
    res[2 * 1:n_x] = y

    res
}

str_trim = function(x, n_first = 0, n_last = 0){
    # positive values: first
    # negative values: last

    if(is.character(n_first)){
        n_first = nchar(n_first)
    } else if(n_first < 0){
        n_last = -n_first
        n_first = 0
    }

    res = substr(x, 1 + n_first, nchar(x) - n_last)
}

str_split = function(x, split){
    # Simple wrapper

    fixed = TRUE
    perl = FALSE
    if(grepl("@", split, fixed = TRUE)){
        if(grepl("^\\\\@", split)){
            split = str_trim(split, 1)
        } else if(grepl("^@", split)){
            split = str_trim(split, 1)
            fixed = FALSE
            perl = TRUE
        }
    }

    strsplit(x, split, fixed = fixed, perl = perl)
}

NA_fun = function(..., df){
    dots = list(...)
    n_dots = length(dots)

    if(n_dots == 0){
        res = !complete.cases(df)
        return(res)
    }

    res = is.na(dots[[1]])

    if(n_dots == 1){
        return(res)
    }

    for(i in 2:n_dots){
        res = res | is.na(dots[[i]])
    }
    res
}


eval_dot = function(x, up = 1){
    # Note that the right use of up is ESSENTIAL
    # if must refer to the parent frame from which the main call has been made
    # Otherwise an error will be thrown of "." not existing

    x_dp = deparse(substitute(x), 300)

    sysOrigin = sys.parent(up)
    mc = match.call(sys.function(sysOrigin), sys.call(sysOrigin))

    if(!x_dp %in% names(mc)){
        return(x)
    }

    my_list = list("." = list)

    eval(mc[[x_dp]], my_list, parent.frame(up + 1))
}


is_user_level_call = function(){
    length(sys.calls()) <= 2
}


is_calling_fun = function(pattern){
    sc_all = sys.calls()
    n_sc = length(sc_all)
    if(n_sc > 2){

        if(grepl(".fixest", sc_all[[n_sc - 1]][[1]], fixed = TRUE)){
            if(n_sc == 3){
                return(FALSE)
            }

            sc = sc_all[[n_sc - 3]]
        } else {
            sc = sc_all[[n_sc - 2]]
        }

        return(grepl(pattern, as.character(sc[[1]])))
    }

    FALSE
}
