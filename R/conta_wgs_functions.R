#' @useDynLib conta
#' @importFrom Rcpp sourceCpp
NULL

#' Clean given object and perform garbage collection
#'
#' @param object to be cleaned
#'
#' @export
clean <- function(object) {
  rm(object)
  gc()
}

#' Return nucleotide bases
#'
#' @export
getBases <- function() {
  return (c("A", "T", "G", "C"))
}
