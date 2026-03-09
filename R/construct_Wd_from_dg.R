#' Construct a quadratic desired-gain matrix `W_d` automatically for QGSI
#'
#' @description
#' Builds a symmetric quadratic desired-gain matrix for the desired-gain QGSI formulation
#' \deqn{QGSI\_DG_i = d' \hat{\gamma}_i + \hat{\gamma}_i' W_d \hat{\gamma}_i.}
#'
#' @param gebv_data A `data.frame` or `data.table` containing trait GEBVs.
#' @param trait_cols Character vector of trait columns.
#' @param dg Named numeric desired-gain vector.
#' @param lower_is_better Optional traits where lower values are preferred.
#' @param center_traits Logical; if `TRUE`, center GEBVs before deriving trait relationships.
#' @param scale_traits Logical; if `TRUE`, scale GEBVs before deriving trait relationships.
#' @param impute_missing Logical; if `TRUE`, impute missing values by trait mean.
#' @param method One of `"hybrid"`, `"diag_only"`, `"dg_outer"`, or `"corr_weighted"`.
#' @param base_diag Optional named baseline diagonal weights.
#' @param lambda_diag Numeric weight for the diagonal component.
#' @param lambda_outer Numeric weight for the desired-gain outer-product component.
#' @param lambda_corr Numeric weight for the correlation-weighted component.
#' @param corr_power Numeric exponent applied to absolute correlations.
#' @param offdiag_shrink Numeric shrinkage factor for off-diagonal terms.
#' @param positive_diagonal_only Logical; if `TRUE`, force non-negative diagonals.
#' @param normalize One of `"trace"`, `"max_abs"`, or `"none"`.
#' @param debug Logical; if `TRUE`, print debug messages.
#'
#' @return A list containing `W_d` and its construction components.
#'
#' @import data.table
#' @export
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