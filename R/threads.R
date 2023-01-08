# TODO: OMP functions
check_set_nthreads <- function(nthreads) {
  # Simple function that checks that the number of threads is valid
  set_up(1)

  check_value(nthreads, "integer scalar GE{0} | numeric scalar GT{0} LT{1}", .message = paste0("The argument 'nthreads' must be an integer lower or equal to the number of threads available (", max(cpp_get_nb_threads(), 1), "). It can be equal to 0 which means all threads. Alternatively, if equal to a number strictly between 0 and 1, it represents the fraction of all threads to be used."))

  max_threads <- cpp_get_nb_threads()

  # # To add later
  # if(cpp_is_in_fork()) return(1)

  if (nthreads == 0) {
    nthreads <- max(max_threads, 1)
  } else if (nthreads >= 1) {
    if (max_threads == 0) {
      warn_up("OpenMP not detected: cannot use ", nthreads, " threads, single-threaded mode instead.")
      nthreads <- 1
    } else if (nthreads > max_threads) {
      warn_up("Asked for ", nthreads, " threads while the maximum is ", max_threads, ". Set to ", max_threads, " threads instead.")
      nthreads <- max_threads
    }
  }

  nthreads
}

#' Sets/gets the number of threads to use in \code{fixest} functions
#'
#' Sets/gets the default number of threads to used in \code{fixest} estimation functions. The default is the maximum number of threads minus two.
#'
#'
#'
#' @param nthreads The number of threads. Can be: a) an integer lower than, or equal to, the maximum number of threads; b) 0: meaning all available threads will be used; c) a number strictly between 0 and 1 which represents the fraction of all threads to use. If missing, the default is to use 50\% of all threads.
#' @param save Either a logical or equal to \code{"reset"}. Default is \code{FALSE}. If \code{TRUE} then the value is set permanently at the project level, this means that if you restart R, you will still obtain the previously saved defaults. This is done by writing in the \code{".Renviron"} file, located in the project's working directory, hence we must have write permission there for this to work, and only works with Rstudio. If equal to "reset", the default at the project level is erased. Since there is writing in a file involved, permission is asked to the user.
#'
#' @author
#' Laurent Berge
#'
#'
#' @examples
#'
#' # Gets the current number of threads
#' getFixest_nthreads()
#'
#' # To set multi-threading off:
#' setFixest_nthreads(1)
#'
#' # To set it back to default:
#' setFixest_nthreads()
#'
setFixest_nthreads <- function(nthreads = 1, save = FALSE) {
  # TODO: OMP functions
  # By default, we use only 50% of threads (never use all)
  max_CRAN <- as.numeric(Sys.getenv("OMP_THREAD_LIMIT"))
  max_CRAN[is.na(max_CRAN)] <- 1000
  max_threads <- min(cpp_get_nb_threads(), 1000, max_CRAN) # we cap at 1k nthreads

  check_arg_plus(save, "logical scalar | match(reset)")

  do_reset <- identical(save, "reset")

  # TODO: OMP functions
  if (missing(nthreads) || is.null(nthreads)) {
    # We first get the default from the environment variable
    # Original behaviour: If it is missing => 50% of all threads
    # nthreads_default = 0.5 50% of all available threads (usually equiv to the nber of procs)

    # New behaviour: detect number of cores first, and then set floor(cores / 2)

    nthreads_default <- renvir_get("fixest_nthreads")

    nthreads_floor <- as.integer(floor(cpp_get_nb_threads() / 2))

    if (!do_reset && !is.null(nthreads_default)) {
      if (!isScalar(nthreads_default) || nthreads_default < 0) {
        warning("The variable setting the number of threads in the .Renviron file is corrupted. It's value has been reset.")
        renvir_update("fixest_nthreads", NULL)
        nthreads_default <- nthreads_floor
      }
    } else {
      nthreads_default <- nthreads_floor
    }

    nthreads <- check_set_nthreads(nthreads_default)
  }

  nthreads <- check_set_nthreads(nthreads)

  if (do_reset) {
    renvir_update("fixest_nthreads", NULL)
  } else if (save) {
    renvir_update("fixest_nthreads", nthreads)
  }

  options("fixest_nthreads" = nthreads)

  invisible()
}

#' @rdname setFixest_nthreads
getFixest_nthreads <- function() {
  x <- getOption("fixest_nthreads")
  if (length(x) != 1 || !is.numeric(x) || is.na(x) || x %% 1 != 0 || x < 0) {
    stop("The value of getOption(\"fixest_nthreads\") is currently not legal. Please use function setFixest_nthreads to set it to an appropriate value. ")
  }

  x
}
