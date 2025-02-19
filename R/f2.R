#' Compute Jaccard Index Using Parallel Rcpp
#'
#' This function takes a matrix, removes rows where all values are zero (or NA),
#' and computes the Jaccard index using an Rcpp parallel function.
#'
#'
#' @importFrom RcppParallel RcppParallelLibs
#' @param smats A numeric matrix where Jaccard Index will be computed between the columns.
#' @return A matrix containing the Jaccard index between columns.
#' @importFrom stats as.dist
#' @export

msca <- function(smats) {
  Rcpp::sourceCpp("./src/FAST_JACK.cpp")
  # Ensure the input is a numeric matrix
  if (!is.matrix(smats) || !is.numeric(smats)) {
    stop("Input must be a numeric matrix.")
  }

  # Remove rows where all values are zero or NA
  i1 <- rowSums(smats, na.rm = TRUE) != 0
  mat2 <- smats[i1, , drop = FALSE]  # Ensure it remains a matrix even if 1 column/row

  # Get patient ids
  pid <- colnames( mat2 )

  # Compute the Jaccard index using the Rcpp parallel function
  d2.1 <- jaccard_index_rcpp_parallel(mat2, mat2)

  dimnames(d2.1) <- list( pid , pid )
  d2 <- as.dist(d2.1)
  # Return the result
  return(d2)
}
