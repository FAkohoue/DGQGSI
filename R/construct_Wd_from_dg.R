#' Construct a quadratic desired-gain matrix `W_d` automatically for QGSI
#'
#' @description
#' Builds a symmetric quadratic desired-gain matrix for the desired-gain QGSI
#' formulation:
#' \deqn{QGSI\_DG_i = d' \hat{\gamma}_i + \hat{\gamma}_i' W_d \hat{\gamma}_i}
#'
#' This helper function provides a flexible way to construct the quadratic matrix
#' `W_d` from a desired-gain vector and the observed trait relationships in a
#' GEBV dataset.
#'
#' @details
#' The matrix `W_d` governs the quadratic part of the desired-gain QGSI. It can
#' be interpreted as controlling how strongly:
#' \itemize{
#'   \item squared trait effects contribute to the index,
#'   \item pairwise trait interactions contribute to the index,
#'   \item correlation structure among traits influences the quadratic term.
#' }
#'
#' The function supports four strategies for constructing `W_d`:
#' \describe{
#'   \item{`"diag_only"`}{Builds a purely diagonal matrix, so the quadratic term
#'   contains only squared trait effects and no cross-trait interactions.}
#'
#'   \item{`"dg_outer"`}{Uses the outer product of the desired-gain vector,
#'   \eqn{d d'}, to construct `W_d`. This introduces pairwise interaction terms
#'   based only on desired gains.}
#'
#'   \item{`"corr_weighted"`}{Constructs `W_d` by weighting the desired-gain
#'   outer product by the empirical trait correlation matrix derived from the
#'   transformed GEBV data. This allows the quadratic matrix to reflect observed
#'   relationships among traits.}
#'
#'   \item{`"hybrid"`}{Combines a diagonal component, an outer-product component,
#'   and a correlation-weighted component. This is the most flexible option and
#'   is generally a good default when both squared effects and trait
#'   interactions are of interest.}
#' }
#'
#' The construction workflow is:
#' \enumerate{
#'   \item convert selected GEBV columns to numeric;
#'   \item optionally impute missing values;
#'   \item optionally flip traits in `lower_is_better`;
#'   \item optionally center and/or scale the GEBV matrix;
#'   \item compute a trait correlation matrix from the transformed GEBV matrix;
#'   \item construct the diagonal, outer-product, and correlation-weighted
#'         components;
#'   \item combine these components according to `method` and the lambda weights;
#'   \item optionally shrink off-diagonal terms;
#'   \item enforce symmetry and optional diagonal positivity;
#'   \item optionally normalize the final matrix.
#' }
#'
#' If trait GEBVs are already centered and scaled, you may wish to use:
#' \itemize{
#'   \item `center_traits = FALSE`
#'   \item `scale_traits = FALSE`
#' }
#'
#' The returned `W_d` can then be passed directly to
#' [run_qgsi_desired_gain()].
#'
#' @param gebv_data A `data.frame` or `data.table` containing trait GEBVs.
#'   Each row should represent one genotype, and the columns listed in
#'   `trait_cols` should contain numeric or coercible-to-numeric trait GEBVs.
#'
#' @param trait_cols Character vector giving the names of the trait columns to
#'   use when constructing `W_d`. These columns must exist in `gebv_data`.
#'   Their order determines the order used in the returned matrix.
#'
#' @param dg Named numeric desired-gain vector, with names matching
#'   `trait_cols`. These desired gains are used to build the outer-product and
#'   correlation-weighted components of `W_d`, and to define default diagonal
#'   weights when `base_diag = NULL`.
#'
#' @param lower_is_better Optional character vector of trait names for which
#'   smaller values are favorable. These traits are internally multiplied by
#'   `-1` so that all traits are oriented in a common favorable direction before
#'   correlations and matrix components are computed.
#'
#' @param center_traits Logical. If `TRUE`, center the selected GEBV columns by
#'   subtracting their means before deriving trait relationships. This is useful
#'   when trait GEBVs are not already centered around zero.
#'
#' @param scale_traits Logical. If `TRUE`, divide each selected trait column by
#'   its standard deviation after optional centering. This is useful when traits
#'   are on different scales and you want the correlation structure to be based
#'   on standardized variables.
#'
#' @param impute_missing Logical. If `TRUE`, missing values are imputed by the
#'   mean of the corresponding trait column. If `FALSE`, the function stops with
#'   an error if missing values are present in the selected trait columns.
#'
#' @param method Character string specifying how `W_d` should be constructed.
#'   One of:
#'   \itemize{
#'     \item `"hybrid"`
#'     \item `"diag_only"`
#'     \item `"dg_outer"`
#'     \item `"corr_weighted"`
#'   }
#'
#' @param base_diag Optional named numeric vector giving baseline diagonal
#'   weights. Names must match `trait_cols`. If `NULL`, the function uses
#'   `abs(dg)` as the default diagonal weights.
#'
#' @param lambda_diag Numeric multiplier applied to the diagonal component of
#'   `W_d`. Larger values increase the contribution of squared trait effects.
#'
#' @param lambda_outer Numeric multiplier applied to the desired-gain outer
#'   product component. Larger values increase the contribution of pairwise trait
#'   interactions implied by the desired-gain vector alone.
#'
#' @param lambda_corr Numeric multiplier applied to the correlation-weighted
#'   component. Larger values increase the influence of the empirical trait
#'   correlation structure on `W_d`.
#'
#' @param corr_power Numeric exponent applied to the absolute trait correlations
#'   before building the correlation-weighted component. Values greater than `1`
#'   emphasize strong correlations more heavily; values between `0` and `1`
#'   reduce contrast among correlation magnitudes.
#'
#' @param offdiag_shrink Numeric shrinkage factor applied to off-diagonal
#'   elements of `W_d`. Must be between `0` and `1`. A value of `1` leaves
#'   off-diagonal elements unchanged; smaller values weaken cross-trait
#'   interaction terms.
#'
#' @param positive_diagonal_only Logical. If `TRUE`, negative diagonal entries
#'   are replaced by zero after construction. This ensures non-negative diagonal
#'   contributions in the final matrix.
#'
#' @param normalize Character string controlling normalization of the final
#'   matrix. One of:
#'   \itemize{
#'     \item `"trace"`: divide `W_d` by the sum of absolute diagonal values;
#'     \item `"max_abs"`: divide `W_d` by its maximum absolute entry;
#'     \item `"none"`: do not normalize.
#'   }
#'
#' @param debug Logical. If `TRUE`, print progress and diagnostic messages during
#'   matrix construction.
#'
#' @return A list containing:
#' \describe{
#'   \item{W_d}{The final symmetric quadratic desired-gain matrix.}
#'   \item{trait_matrix_used}{A `data.table` containing the transformed GEBV
#'     matrix actually used to derive trait relationships after optional
#'     imputation, direction flipping, centering, and scaling.}
#'   \item{cor_matrix}{The empirical trait correlation matrix computed from the
#'     transformed GEBV matrix.}
#'   \item{diag_component}{The diagonal component used in constructing `W_d`.}
#'   \item{outer_component}{The desired-gain outer-product component.}
#'   \item{corr_component}{The correlation-weighted component.}
#'   \item{method}{The construction method used.}
#' }
#'
#' @examples
#' trait_cols <- c("YLD", "MY", "MI", "BL", "NBL", "VHB")
#'
#' dg <- c(
#'   YLD = 1.5,
#'   MY  = 0.5,
#'   MI  = 0.5,
#'   BL  = 1.0,
#'   NBL = 1.0,
#'   VHB = 1.0
#' )
#'
#' ext <- system.file("extdata", package = "DGQGSI")
#' gebv <- data.table::fread(file.path(ext, "example_gebv.csv"))
#'
#' Wobj <- construct_Wd_from_dg(
#'   gebv_data = gebv[, c("GenoID", trait_cols), with = FALSE],
#'   trait_cols = trait_cols,
#'   dg = dg,
#'   lower_is_better = c("BL", "NBL", "VHB"),
#'   center_traits = FALSE,
#'   scale_traits = FALSE,
#'   method = "hybrid",
#'   lambda_diag = 1.0,
#'   lambda_outer = 0.5,
#'   lambda_corr = 1.0,
#'   corr_power = 1,
#'   offdiag_shrink = 0.5,
#'   normalize = "trace",
#'   debug = FALSE
#' )
#'
#' Wobj$W_d
#' Wobj$cor_matrix
#'
#' @references
#' Cerón-Rojas JJ, Montesinos-López OA, Montesinos-López A, et al. (2026).
#' Nonlinear genomic selection index accelerates multi-trait crop improvement.
#' \emph{Nature Communications}, 17:1991.
#' doi:10.1038/s41467-026-69890-3
#'
#' @import data.table
#' @export
#' 
construct_Wd_from_dg <- function(
    gebv_data,
    trait_cols,
    dg,
    lower_is_better = NULL,
    center_traits = TRUE,
    scale_traits = TRUE,
    impute_missing = TRUE,
    method = c("hybrid", "diag_only", "dg_outer", "corr_weighted"),
    base_diag = NULL,
    lambda_diag = 1.0,
    lambda_outer = 0.5,
    lambda_corr = 1.0,
    corr_power = 1,
    offdiag_shrink = 1.0,
    positive_diagonal_only = TRUE,
    normalize = c("trace", "max_abs", "none"),
    debug = TRUE
) {
  method <- match.arg(method)
  normalize <- match.arg(normalize)
  
  .dgqgsi_dbg(debug, "============================================================")
  .dgqgsi_dbg(debug, "Starting construct_Wd_from_dg()")
  .dgqgsi_dbg(debug, "Method: %s", method)
  
  dt <- data.table::as.data.table(data.table::copy(gebv_data))
  
  if (!all(trait_cols %in% names(dt))) {
    miss <- setdiff(trait_cols, names(dt))
    stop("Missing trait columns in gebv_data: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  
  dg <- dg[trait_cols]
  if (anyNA(dg)) {
    stop("dg missing for traits: ", paste(trait_cols[is.na(dg)], collapse = ", "), call. = FALSE)
  }
  dg <- as.numeric(dg)
  names(dg) <- trait_cols
  
  for (tr in trait_cols) {
    dt[, (tr) := as.numeric(get(tr))]
  }
  
  if (isTRUE(impute_missing)) {
    for (tr in trait_cols) {
      if (anyNA(dt[[tr]])) {
        mu <- mean(dt[[tr]], na.rm = TRUE)
        dt[is.na(get(tr)), (tr) := mu]
      }
    }
  } else if (anyNA(as.matrix(dt[, trait_cols, with = FALSE]))) {
    stop("Missing values found in gebv_data while impute_missing = FALSE.", call. = FALSE)
  }
  
  if (!is.null(lower_is_better)) {
    if (!all(lower_is_better %in% trait_cols)) {
      bad <- setdiff(lower_is_better, trait_cols)
      stop("These lower_is_better traits are not in trait_cols: ", paste(bad, collapse = ", "), call. = FALSE)
    }
    dt[, (lower_is_better) := lapply(.SD, function(x) -x), .SDcols = lower_is_better]
  }
  
  X <- as.matrix(dt[, trait_cols, with = FALSE])
  
  if (isTRUE(center_traits)) {
    X <- sweep(X, 2, colMeans(X), FUN = "-")
  }
  
  if (isTRUE(scale_traits)) {
    sds <- apply(X, 2, stats::sd)
    bad_sd <- !is.finite(sds) | sds == 0
    sds[bad_sd] <- 1
    X <- sweep(X, 2, sds, FUN = "/")
  }
  
  colnames(X) <- trait_cols
  
  R <- stats::cor(X, use = "pairwise.complete.obs")
  R <- 0.5 * (R + t(R))
  diag(R) <- 1
  
  if (is.null(base_diag)) {
    diag_vec <- abs(dg)
  } else {
    base_diag <- base_diag[trait_cols]
    if (anyNA(base_diag)) {
      stop("base_diag missing for traits: ", paste(trait_cols[is.na(base_diag)], collapse = ", "), call. = FALSE)
    }
    diag_vec <- as.numeric(base_diag)
    names(diag_vec) <- trait_cols
  }
  
  diag_component <- diag(diag_vec, nrow = length(trait_cols), ncol = length(trait_cols))
  rownames(diag_component) <- colnames(diag_component) <- trait_cols
  
  outer_component <- tcrossprod(dg)
  rownames(outer_component) <- colnames(outer_component) <- trait_cols
  
  corr_weight_matrix <- sign(R) * (abs(R)^corr_power)
  corr_component <- outer_component * corr_weight_matrix
  rownames(corr_component) <- colnames(corr_component) <- trait_cols
  
  if (method == "diag_only") W_d <- lambda_diag * diag_component
  if (method == "dg_outer") W_d <- lambda_outer * outer_component
  if (method == "corr_weighted") W_d <- lambda_corr * corr_component
  if (method == "hybrid") {
    W_d <- lambda_diag * diag_component +
      lambda_outer * outer_component +
      lambda_corr * corr_component
  }
  
  if (!is.numeric(offdiag_shrink) || length(offdiag_shrink) != 1 || offdiag_shrink < 0 || offdiag_shrink > 1) {
    stop("offdiag_shrink must be a single number between 0 and 1.", call. = FALSE)
  }
  
  if (offdiag_shrink < 1) {
    offdiag_mask <- row(W_d) != col(W_d)
    W_d[offdiag_mask] <- W_d[offdiag_mask] * offdiag_shrink
  }
  
  W_d <- 0.5 * (W_d + t(W_d))
  
  if (isTRUE(positive_diagonal_only)) {
    diag(W_d) <- pmax(diag(W_d), 0)
  }
  
  if (normalize == "max_abs") {
    mx <- max(abs(W_d))
    if (is.finite(mx) && mx > 0) W_d <- W_d / mx
  }
  
  if (normalize == "trace") {
    tr <- sum(diag(abs(W_d)))
    if (is.finite(tr) && tr > 0) W_d <- W_d / tr
  }
  
  list(
    W_d = W_d,
    trait_matrix_used = data.table::as.data.table(X),
    cor_matrix = R,
    diag_component = diag_component,
    outer_component = outer_component,
    corr_component = corr_component,
    method = method
  )
}