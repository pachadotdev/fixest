parse_int <- function(x) {
  if (!is.integer(x)) x <- as.integer(x)
  x
}

parse_chr <- function(x) {
  if (!is.character(x)) x <- as.character(x)
  x
}

parse_int_mat <- function(x) {
  if (!(storage.mode(x) == "integer")) storage.mode(x) <- "integer"
  x
}
