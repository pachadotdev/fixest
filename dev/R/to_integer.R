#' Fast transform of any type of vector(s) into an integer vector
#'
#' Tool to transform any type of vector, or even combination of vectors, into an integer vector ranging from 1 to the number of unique values. This actually creates an unique identifier vector.
#'
#' @param ... Vectors of any type, to be transformed in integer.
#' @param sorted Logical, default is `FALSE`. Whether the integer vector should make reference to sorted values?
#' @param add_items Logical, default is `FALSE`. Whether to add the unique values of the original vector(s). If requested, an attribute `items` is created containing the values (alternatively, they can appear in a list if `items.list=TRUE`).
#' @param items.list Logical, default is `FALSE`. Only used if `add_items=TRUE`. If `TRUE`, then a list of length 2 is returned with `x` the integer vector and `items` the vector of items.
#' @param multi.df Logical, default is `FALSE`. If `TRUE` then a data.frame listing the unique elements is returned in the form of a data.frame. Ignored if `add_items = FALSE`.
#' @param multi.join Character scalar used to join the items of multiple vectors. The default is `"_"`. Ignored if `add_items = FALSE`.
#' @param internal Logical, default is `FALSE`. For programming only. If this function is used within another function, setting `internal = TRUE` is needed to make the evaluation of `...` valid. End users of `to_integer` should not care.
#'
#'
#' @return
#' Reruns a vector of the same length as the input vectors.
#' If `add_items=TRUE` and `items.list=TRUE`, a list of two elements is returned: `x` being the integer vector and `items` being the unique values to which the values in `x` make reference.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' x1 <- iris$Species
#' x2 <- as.integer(iris$Sepal.Length)
#'
#' # transforms the species vector into integers
#' to_integer(x1)
#'
#' # To obtain the "items":
#' to_integer(x1, add_items = TRUE)
#' # same but in list form
#' to_integer(x1, add_items = TRUE, items.list = TRUE)
#'
#' # transforms x2 into an integer vector from 1 to 4
#' to_integer(x2, add_items = TRUE)
#'
#' # To have the sorted items:
#' to_integer(x2, add_items = TRUE, sorted = TRUE)
#'
#' # The result can safely be used as an index
#' res <- to_integer(x2, add_items = TRUE, sorted = TRUE, items.list = TRUE)
#' all(res$items[res$x] == x2)
#'
#' # Multiple vectors ----
#'
#' to_integer(x1, x2, add_items = TRUE)
#'
#' # You can use multi.join to handle the join of the items:
#' to_integer(x1, x2, add_items = TRUE, multi.join = "; ")
#'
#' @export
to_integer <- function(..., sorted = FALSE, add_items = FALSE, items.list = FALSE,
                       multi.df = FALSE, multi.join = "_", internal = FALSE) {
  if (!internal) check_arg(..., "vector mbt")
  check_arg(sorted, add_items, items.list, "logical scalar")
  check_arg(multi.join, "character scalar")

  dots <- list(...)

  # Removing NAs
  Q <- length(dots)
  n_all <- lengths(dots)
  n <- n_all[1]

  if (length(unique(n_all)) != 1) stop("All elements in `...` should be of the same length (current lenghts are ", enumerate_items(n_all), ").")

  is_na <- is.na(dots[[1]])
  for (q in seq(from = 2, length.out = Q - 1)) {
    is_na <- is_na | is.na(dots[[q]])
  }

  ANY_NA <- FALSE
  if (any(is_na)) {
    ANY_NA <- TRUE

    if (all(is_na)) {
      message("NOTE: All values are NA.")
      res <- rep(NA, n)
      if (add_items) {
        if (items.list) {
          res <- list(x = res, items = NA)
        } else {
          attr(res, "items") <- NA
        }
      }

      return(res)
    }

    for (q in 1:Q) dots[[q]] <- dots[[q]][!is_na]
  }

  #
  # Creating the ID
  #

  if (Q == 1) {
    if (sorted && !is.numeric(dots[[1]]) && !is.character(dots[[1]])) {
      # general way => works for any type with a sort method
      f <- dots[[1]]
      res_raw <- quickUnclassFactor(f, addItem = TRUE, sorted = FALSE)
      obs_1st <- cpp_get_first_item(res_raw$x, length(res_raw$items))
      f_unik <- f[obs_1st]
      f_order <- order(f_unik)
      x_new <- order(f_order)[res_raw$x]
      if (add_items) {
        items_new <- f_unik[f_order]
        res <- list(x = x_new, items = items_new)
      } else {
        res <- x_new
      }
    } else {
      res <- quickUnclassFactor(dots[[1]], addItem = add_items, sorted = sorted)
    }
  } else {
    QUF_raw <- list()
    for (q in 1:Q) {
      QUF_raw[[q]] <- quickUnclassFactor(dots[[q]], sorted = FALSE, addItem = TRUE)
    }

    # Then we combine
    power <- floor(1 + log10(sapply(QUF_raw, function(x) length(x$items))))

    is_large <- sum(power) > 14
    if (is_large) {
      # 15 Aug 2021, finally found a solution. It was so obvious with hindsight...
      QUF_raw_value <- lapply(QUF_raw, `[[`, 1)
      order_index <- do.call(order, QUF_raw_value)
      index <- cpp_combine_clusters(QUF_raw_value, order_index)
    } else {
      # quicker, but limited by the precision of doubles
      index <- QUF_raw[[1]]$x
      for (q in 2:Q) {
        index <- index + QUF_raw[[q]]$x * 10**sum(power[1:(q - 1)])
      }
    }

    res <- quickUnclassFactor(index, addItem = add_items || sorted, sorted = sorted)

    if (add_items || sorted) {
      # we re order appropriately
      # f prefix means factor

      obs_1st <- cpp_get_first_item(res$x, length(res$items))

      f_all <- list()
      for (q in 1:Q) {
        f_all[[q]] <- dots[[q]][obs_1st]
      }

      f_order <- do.call("order", f_all)

      x_new <- order(f_order)[res$x]

      if (multi.df) {
        # Putting into a DF => we take care of names
        mc_dots <- match.call(expand.dots = FALSE)[["..."]]
        n_dots <- length(mc_dots)
        mc_dots_names <- names(mc_dots)
        if (is.null(mc_dots_names)) mc_dots_names <- character(n_dots)

        my_names <- character(n_dots)
        for (q in 1:n_dots) {
          if (nchar(mc_dots_names[q]) > 0) {
            my_names[q] <- mc_dots_names[q]
          } else {
            my_names[q] <- deparse_long(mc_dots[[q]])
          }
        }

        names(f_all) <- my_names

        f_df <- as.data.frame(f_all)
        items_new <- f_df[f_order, , drop = FALSE]
        row.names(items_new) <- seq_len(nrow(items_new))
      } else {
        # we "paste" them
        arg_list <- f_all
        arg_list$sep <- multi.join
        f_char <- do.call("paste", arg_list)
        items_new <- f_char[f_order]
      }

      if (add_items) {
        res <- list(x = x_new, items = items_new)
      } else {
        res <- x_new
      }
    }
  }

  if (ANY_NA) {
    if (is.list(res)) {
      x <- res$x
    } else {
      x <- res
    }

    x_na <- rep(NA, n)
    x_na[!is_na] <- x

    if (is.list(res)) {
      res$x <- x_na
    } else {
      res <- x_na
    }
  }

  if (add_items && isFALSE(items.list)) {
    res_tmp <- res$x
    attr(res_tmp, "items") <- res$items
    res <- res_tmp
  }

  res
}
