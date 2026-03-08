# Internal utilities for DGQGSI

# nocov start
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "..dg_cols_present",
    "..qgsi_cols_present",
    "SelectionIndex",
    "Selected",
    "Rank_SelectionIndex",
    "LinearDGPart",
    "QuadraticDGPart",
    "QGSI_DG",
    "Rank_QGSI_DG",
    "Rank_LinearDGPart",
    "Rank_QuadraticDGPart",
    "DG_SelectionIndex",
    "DG_Rank_SelectionIndex",
    "DG_Selected",
    "RankDiff_DG_minus_QGSI",
    "AbsRankDiff_DG_vs_QGSI",
    "SignSame_DG_QGSI"
  ))
}
# nocov end

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
