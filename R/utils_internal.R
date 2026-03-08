# Internal utilities for DGQGSI

#' Internal debug message helper
#'
#' @param debug Logical.
#' @param ... Arguments passed to `sprintf()`.
#'
#' @keywords internal
.dgqgsi_dbg <- function(debug, ...) {
  if (isTRUE(debug)) message(sprintf(...))
}

#' Internal z-score helper
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector.
#' @keywords internal
.dgqgsi_z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

#' Validate a square matrix
#'
#' @param M Matrix.
#' @param p Expected dimension.
#' @param nm Object name for messages.
#'
#' @keywords internal
.validate_square_matrix <- function(M, p, nm = "matrix") {
  if (!is.matrix(M)) stop(sprintf("%s must be a matrix.", nm), call. = FALSE)
  if (nrow(M) != p || ncol(M) != p) {
    stop(sprintf("%s must be a %d x %d matrix.", nm, p, p), call. = FALSE)
  }
  invisible(TRUE)
}
