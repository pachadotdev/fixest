# setters/getters -----

#' Sets/gets whether to display notes in \code{fixest} estimation functions
#'
#' Sets/gets the default values of whether notes (informing for NA and observations removed) should be displayed in \code{fixest} estimation functions.
#'
#' @param x A logical. If \code{FALSE}, then notes are permanently removed.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Change default with
#' setFixest_notes(FALSE)
#' feols(Ozone ~ Solar.R, airquality)
#'
#' # Back to default which is TRUE
#' setFixest_notes(TRUE)
#' feols(Ozone ~ Solar.R, airquality)
#'
setFixest_notes = function(x){
    check_arg(x, "mbt logical scalar")

    options("fixest_notes" = x)
}

#' @rdname setFixest_notes
getFixest_notes = function(){

    x = getOption("fixest_notes")
    if(length(x) != 1 || !is.logical(x) || is.na(x)){
        stop("The value of getOption(\"fixest_notes\") is currently not legal. Please use function setFixest_notes to set it to an appropriate value. ")
    }

    x
}

#' Sets/gets the dictionary relabeling the variables
#'
#' Sets/gets the default dictionary used in the function \code{\link[fixest]{etable}}, \code{\link[fixest]{did_means}} and \code{\link[fixest]{coefplot}}. The dictionaries are used to relabel variables (usually towards a fancier, more explicit formatting) when exporting them into a Latex table or displaying in graphs. By setting the dictionary with \code{setFixest_dict}, you can avoid providing the argument \code{dict}.
#'
#'
#' @param dict A named character vector. E.g. to change my variable named "a" and "b" to (resp.) "$log(a)$" and "$bonus^3$", then use \code{dict = c(a="$log(a)$", b3="$bonus^3$")}. This dictionary is used in Latex tables or in graphs by the function \code{\link[fixest]{coefplot}}. If you want to separate Latex rendering from rendering in graphs, use an ampersand first to make the variable specific to \code{coefplot}.
#'
#' @author
#' Laurent Berge
#'
#'
#' @examples
#'
#' data(trade)
#' est = feols(log(Euros) ~ log(dist_km)|Origin+Destination+Product, trade)
#' # we export the result & rename some variables
#' esttex(est, dict = c("log(Euros)"="Euros (ln)", Origin="Country of Origin"))
#'
#' # If you export many tables, it can be more convenient to use setFixest_dict:
#' setFixest_dict(c("log(Euros)"="Euros (ln)", Origin="Country of Origin"))
#' esttex(est) # variables are properly relabeled
#'
setFixest_dict = function(dict){

    if(missing(dict) || is.null(dict)){
        options("fixest_dict" = NULL)
        return(invisible())
    }

    #
    # Controls
    #

    if(!is.character(dict) || !isVector(dict)){
        stop("Argument 'dict' must be a character vector.")
    }

    if(anyNA(dict)){
        stop("Argument 'dict' must be a character vector without NAs.")
    }

    # Formatting the names
    dict_names = names(dict)
    if(is.null(dict_names)){
        stop("Argument 'dict', the dictionary, must be a named vector. Currently it has no names.")
    }

    dict_names = gsub(" +", "", dict_names)
    td = table(dict_names)
    if(any(td > 1)){
        qui = which(dict_names %in% names(td)[td > 1])
        name_dup = unique(names(dict)[qui])
        stop("Argument 'dict' contains duplicated names: ", enumerate_items(name_dup))
    }

    options("fixest_dict" = dict)
}

#' @rdname setFixest_dict
getFixest_dict = function(){

    x = getOption("fixest_dict")
    if(length(x) > 0){
        if(!is.character(x) || !isVector(x) || anyNA(x)){
            stop("The value of getOption(\"fixest_dict\") is currently not legal. Please use function setFixest_dict to set it to an appropriate value. ")
        }
    }

    x
}

#' Sets/gets formula macros
#'
#' You can set formula macros globally with \code{setFixest_fml}. These macros can then be used in \code{fixest} estimations or when using the function \code{\link[fixest:setFixest_fml]{xpd}}.
#'
#' @inherit xpd examples
#'
#' @param ... Definition of the macro variables. Each argument name corresponds to the name of the macro variable. It is required that each macro variable name starts with two dots (e.g. \code{..ctrl}). The value of each argument must be a one-sided formula or a character vector, it is the definition of the macro variable. Example of a valid call: \code{setFixest_fml(..ctrl = ~ var1 + var2)}. In the function \code{xpd}, the default macro variables are taken from \code{getFixest_fml}, any variable in \code{...} will replace these values. You can enclose values in \code{.[]}, if so they will be evaluated from the current environment. For example \code{..ctrl = ~ x.[1:2] + .[z]} will lead to \code{~x1 + x2 + var} if \code{z} is equal to \code{"var"}.
#' @param reset A logical scalar, defaults to \code{FALSE}. If \code{TRUE}, all macro variables are first reset (i.e. deleted).
#'
#' @details
#' In \code{xpd}, the default macro variables are taken from \code{getFixest_fml}. Any value in the \code{...} argument of \code{xpd} will replace these default values.
#'
#' The definitions of the macro variables will replace in verbatim the macro variables. Therefore, you can include multipart formulas if you wish but then beware of the order the the macros variable in the formula. For example, using the airquality data, say you want to set as controls the variable \code{Temp} and \code{Day} fixed-effects, you can do \code{setFixest_fml(..ctrl = ~Temp | Day)}, but then \code{feols(Ozone ~ Wind + ..ctrl, airquality)} will be quite different from \code{feols(Ozone ~ ..ctrl + Wind, airquality)}, so beware!
#'
#' @return
#' The function \code{getFixest_fml()} returns a list of character strings, the names corresponding to the macro variable names, the character strings corresponding to their definition.
#'
#' @seealso
#' \code{\link[fixest]{xpd}} to make use of formula macros.
#'
#'
#'
setFixest_fml = function(..., reset = FALSE){

    check_arg(reset, get, "logical scalar")

    fml_macro = parse_macros(..., reset = reset, frame = parent.frame())

    options("fixest_fml_macro" = fml_macro)

}

#' @rdname setFixest_fml
getFixest_fml = function(){
    fml_macro = getOption("fixest_fml_macro")
    if(is.null(fml_macro)){
        options("fixest_fml_macro" = list())
        fml_macro = list()
    } else if(!is.list(fml_macro)){
        options("fixest_fml_macro" = list())
        warning("The value of getOption(\"fixest_fml_macro\") is not legal, it has been reset. Please use only 'setFixest_fml' to set formula macro variables.")
        fml_macro = list()
    }

    fml_macro
}



#' Default arguments for fixest estimations
#'
#' This function sets globally the default arguments of fixest estimations.
#'
#' @inheritParams feols
#' @inheritParams feNmlm
#' @inheritParams feglm
#'
#' @param reset Logical, default to \code{FALSE}. Whether to reset all values.
#'
#' @return
#' The function \code{getFixest_estimation} returns the currently set global defaults.
#'
#' @examples
#'
#' #
#' # Example: removing singletons is FALSE by default
#' #
#'
#' # => changing this default
#'
#' # Let's create data with singletons
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#' base$fe_singletons = as.character(base$species)
#' base$fe_singletons[1:5] = letters[1:5]
#'
#' res          = feols(y ~ x1 + x2 | fe_singletons, base)
#' res_noSingle = feols(y ~ x1 + x2 | fe_singletons, base, fixef.rm = "single")
#'
#' # New defaults
#' setFixest_estimation(fixef.rm = "single")
#' res_newDefault = feols(y ~ x1 + x2 | fe_singletons, base)
#'
#' etable(res, res_noSingle, res_newDefault)
#'
#' # Resetting the defaults
#' setFixest_estimation(reset = TRUE)
#'
#'
#'
setFixest_estimation = function(data = NULL, panel.id = NULL, fixef.rm = "perfect", fixef.tol = 1e-6, fixef.iter = 10000, collin.tol = 1e-10, lean = FALSE, verbose = 0, warn = TRUE, combine.quick = NULL, demeaned = FALSE, mem.clean = FALSE, glm.iter = 25, glm.tol = 1e-8, reset = FALSE){

    check_arg_plus(fixef.rm, "match(none, perfect, singleton, both)")
    check_arg(fixef.tol, collin.tol, glm.tol, "numeric scalar GT{0}")
    check_arg(fixef.iter, glm.iter, "integer scalar GE{1}")
    check_arg(verbose, "integer scalar GE{0}")
    check_arg(lean, warn, demeaned, mem.clean, reset, "logical scalar")
    check_arg(combine.quick, "NULL logical scalar")
    check_arg(panel.id, "NULL character vector len(,2) no na | os formula")

    if(!missing(data)){
        check_arg(data, "NULL data.frame | matrix")
        if(!is.null(data)){
            data = deparse_long(substitute(data))
            class(data) = "default_data"
        }
    }

    # Getting the existing defaults
    opts = getOption("fixest_estimation")

    if(reset || is.null(opts)){
        opts = list()
    } else if(!is.list(opts)){
        warning("Wrong formatting of option 'fixest_estimation', all options are reset.")
        opts = list()
    }

    # Saving the default values
    mc = match.call()
    args_default = setdiff(names(mc)[-1], "reset")

    # NOTA: we don't allow delayed evaluation => all arguments must have hard values
    for(v in args_default){
        opts[[v]] = eval(as.name(v))
    }

    options(fixest_estimation = opts)

}

#' @rdname setFixest_estimation
getFixest_estimation = function(){
    getOption("fixest_estimation")
}

