#' Desired-gain QGSI using trait GEBVs
#'
#' @description
#' Builds a genomic quadratic selection index for all genotypes using the QGSI
#' structure from the paper, but replacing the linear economic-weight vector
#' with a desired-gain vector:
#' \deqn{QGSI\_DG_i = d' \hat{\gamma}_i + \hat{\gamma}_i' W_d \hat{\gamma}_i.}
#'
#' @details
#' Intended for genomic selection workflows where trait GEBVs are already available.
#' If traits are already standardized and `dg` is in SD units, the recommended
#' settings are `center_traits = FALSE` and `scale_traits = FALSE`.
#'
#' @param init_data A `data.frame` or `data.table` containing genotype IDs and metadata.
#' @param gebv_data A `data.frame` or `data.table` containing genotype IDs and trait GEBVs.
#' @param trait_cols Character vector of trait GEBV column names.
#' @param id_col Character genotype identifier column. Default is `"GenoID"`.
#' @param dg Named numeric desired-gain vector.
#' @param W_d Optional full symmetric quadratic desired-gain matrix.
#' @param quadratic_diag_weights Optional diagonal weights used only when `W_d = NULL`.
#' @param lower_is_better Optional traits where lower values are preferred.
#' @param center_traits Logical; if `TRUE`, center GEBVs before computing the index.
#' @param scale_traits Logical; if `TRUE`, scale GEBVs before computing the index.
#' @param impute_missing Logical; if `TRUE`, impute missing trait values by trait mean.
#' @param return_components Logical; if `TRUE`, return linear and quadratic components.
#' @param debug Logical; if `TRUE`, print debug messages.
#'
#' @return A list containing `W_d`, transformed trait space, ranked genotypes, and component summaries.
#'
#' @import data.table
#' @export
run_qgsi_desired_gain <- function(
    init_data,
    gebv_data,
    trait_cols,
    id_col = "GenoID",
    dg,
    W_d = NULL,
    quadratic_diag_weights = NULL,
    lower_is_better = NULL,
    center_traits = TRUE,
    scale_traits = TRUE,
    impute_missing = TRUE,
    return_components = TRUE,
    debug = TRUE
) {
  dt_init <- data.table::as.data.table(data.table::copy(init_data))
  dt_gebv <- data.table::as.data.table(data.table::copy(gebv_data))
  if (!id_col %in% names(dt_init)) stop(sprintf("id_col '%s' not found in init_data.", id_col), call. = FALSE)
  if (!id_col %in% names(dt_gebv)) stop(sprintf("id_col '%s' not found in gebv_data.", id_col), call. = FALSE)
  if (!all(trait_cols %in% names(dt_gebv))) {
    miss <- setdiff(trait_cols, names(dt_gebv))
    stop("Missing trait columns in gebv_data: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  common_ids <- intersect(dt_init[[id_col]], dt_gebv[[id_col]])
  if (length(common_ids) == 0L) stop("No common IDs between init_data and gebv_data.", call. = FALSE)
  dt_init <- dt_init[get(id_col) %in% common_ids]
  dt_gebv <- dt_gebv[get(id_col) %in% common_ids]
  data.table::setkeyv(dt_init, id_col)
  data.table::setkeyv(dt_gebv, id_col)
  dt_init <- dt_init[dt_gebv[[id_col]]]
  if (!identical(dt_init[[id_col]], dt_gebv[[id_col]])) stop("Alignment failed between init_data and gebv_data.", call. = FALSE)

  for (tr in trait_cols) dt_gebv[, (tr) := as.numeric(get(tr))]
  if (isTRUE(impute_missing)) {
    for (tr in trait_cols) if (anyNA(dt_gebv[[tr]])) dt_gebv[is.na(get(tr)), (tr) := mean(dt_gebv[[tr]], na.rm = TRUE)]
  } else if (anyNA(as.matrix(dt_gebv[, ..trait_cols]))) {
    stop("Missing values detected in gebv_data while impute_missing = FALSE.", call. = FALSE)
  }

  if (!is.null(lower_is_better)) {
    if (!all(lower_is_better %in% trait_cols)) {
      bad <- setdiff(lower_is_better, trait_cols)
      stop("These lower_is_better traits are not in trait_cols: ", paste(bad, collapse = ", "), call. = FALSE)
    }
    dt_gebv[, (lower_is_better) := lapply(.SD, function(x) -x), .SDcols = lower_is_better]
  }

  gebv_mat <- as.matrix(dt_gebv[, ..trait_cols])
  if (isTRUE(center_traits)) gebv_mat <- sweep(gebv_mat, 2, colMeans(gebv_mat), FUN = "-")
  if (isTRUE(scale_traits)) {
    sds <- apply(gebv_mat, 2, stats::sd)
    bad_sd <- !is.finite(sds) | sds == 0
    sds[bad_sd] <- 1
    gebv_mat <- sweep(gebv_mat, 2, sds, FUN = "/")
  }
  colnames(gebv_mat) <- trait_cols

  dg <- dg[trait_cols]
  if (anyNA(dg)) stop("dg missing for traits: ", paste(trait_cols[is.na(dg)], collapse = ", "), call. = FALSE)
  dg <- as.numeric(dg)
  names(dg) <- trait_cols

  p <- length(trait_cols)
  if (is.null(W_d)) {
    if (is.null(quadratic_diag_weights)) {
      quadratic_diag_weights <- rep(1, p)
      names(quadratic_diag_weights) <- trait_cols
    } else {
      quadratic_diag_weights <- quadratic_diag_weights[trait_cols]
      if (anyNA(quadratic_diag_weights)) stop("quadratic_diag_weights missing for traits.", call. = FALSE)
      quadratic_diag_weights <- as.numeric(quadratic_diag_weights)
      names(quadratic_diag_weights) <- trait_cols
    }
    W_d <- diag(quadratic_diag_weights, nrow = p, ncol = p)
    rownames(W_d) <- colnames(W_d) <- trait_cols
  } else {
    .validate_square_matrix(W_d, p, "W_d")
    W_d <- 0.5 * (W_d + t(W_d))
    rownames(W_d) <- colnames(W_d) <- trait_cols
  }

  linear_part <- as.numeric(gebv_mat %*% dg)
  quadratic_part <- apply(gebv_mat, 1, function(x) {
    x <- matrix(x, ncol = 1)
    as.numeric(t(x) %*% W_d %*% x)
  })
  qgsi_dg <- linear_part + quadratic_part

  dt_out <- data.table::copy(dt_init)
  if (isTRUE(return_components)) {
    dt_out[, LinearDGPart := linear_part]
    dt_out[, QuadraticDGPart := quadratic_part]
  }
  dt_out[, QGSI_DG := qgsi_dg]
  dt_out[, Rank_QGSI_DG := data.table::frank(-QGSI_DG, ties.method = "average")]
  if (isTRUE(return_components)) {
    dt_out[, Rank_LinearDGPart := data.table::frank(-LinearDGPart, ties.method = "average")]
    dt_out[, Rank_QuadraticDGPart := data.table::frank(-QuadraticDGPart, ties.method = "average")]
  }
  data.table::setorder(dt_out, -QGSI_DG)

  if (isTRUE(return_components)) {
    component_summary <- data.table::data.table(
      Component = c("LinearDGPart", "QuadraticDGPart", "QGSI_DG"),
      Mean = c(mean(linear_part), mean(quadratic_part), mean(qgsi_dg)),
      SD = c(stats::sd(linear_part), stats::sd(quadratic_part), stats::sd(qgsi_dg)),
      Min = c(min(linear_part), min(quadratic_part), min(qgsi_dg)),
      Max = c(max(linear_part), max(quadratic_part), max(qgsi_dg))
    )
  } else {
    component_summary <- data.table::data.table(
      Component = "QGSI_DG",
      Mean = mean(qgsi_dg), SD = stats::sd(qgsi_dg), Min = min(qgsi_dg), Max = max(qgsi_dg)
    )
  }

  list(
    dg = dg,
    W_d = W_d,
    trait_space_used = data.table::as.data.table(gebv_mat),
    ranked_geno = dt_out,
    component_summary = component_summary
  )
}
