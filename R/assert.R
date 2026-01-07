#' Assert whether an object inherits a certain class
#' 
#' @param x object to be tested
#' @param y class to be tested
#' @note copied from https://github.com/ropensci/taxizedb/blob/master/R/zzz.R
#' @noRd
assert <- function(x, y) {
  if (!is.null(x)) {
    if (!inherits(x, y)) {
      stop(deparse(substitute(x)), " must be of class ",
        paste0(y, collapse = ", "), call. = FALSE)
    }
  }
}