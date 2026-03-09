#' Desired-gain QGSI using trait GEBVs
#'
#' @description
#' Builds a genomic quadratic selection index for all genotypes using the QGSI
#' structure proposed for nonlinear genomic selection, but replacing the linear
#' economic-weight vector with a breeder-defined desired-gain vector.
#'
#' The implemented index is:
#' \deqn{QGSI\_DG_i = d' \hat{\gamma}_i + \hat{\gamma}_i' W_d \hat{\gamma}_i}
#'
#' where:
#' \itemize{
#'   \item \eqn{d} is the desired-gain vector,
#'   \item \eqn{\hat{\gamma}_i} is the vector of trait GEBVs for genotype \eqn{i},
#'   \item \eqn{W_d} is a quadratic desired-gain matrix controlling trait-by-trait
#'         interaction contributions.
#' }
#'
#' @details
#' This function is intended for genomic selection workflows where trait GEBVs
#' are already available for each genotype.
#'
#' The function:
#' \enumerate{
#'   \item aligns `init_data` and `gebv_data` by `id_col`;
#'   \item converts all requested trait columns to numeric;
#'   \item optionally imputes missing values by trait means;
#'   \item optionally flips traits listed in `lower_is_better` so that all traits
#'         are expressed in a common favorable direction;
#'   \item optionally centers and/or scales the trait GEBVs;
#'   \item computes a linear desired-gain component;
#'   \item computes a quadratic interaction component using `W_d`;
#'   \item combines both parts into a final desired-gain QGSI score;
#'   \item ranks all genotypes by the resulting index.
#' }
#'
#' If traits are already standardized and `dg` is already expressed in SD units,
#' the recommended settings are:
#' \itemize{
#'   \item `center_traits = FALSE`
#'   \item `scale_traits = FALSE`
#' }
#'
#' because the trait GEBVs and desired gains are already on the same scale.
#'
#' If `W_d` is not supplied, the function builds a diagonal matrix from
#' `quadratic_diag_weights`. In that case, the quadratic contribution captures
#' only squared trait effects, without off-diagonal trait interactions.
#'
#' @param init_data A `data.frame` or `data.table` containing genotype
#'   identifiers and any metadata to preserve in the final output. Typical
#'   columns may include `GenoID`, family, population, or other descriptors.
#'   The returned table appends QGSI results to this object.
#'
#' @param gebv_data A `data.frame` or `data.table` containing genotype IDs and
#'   trait GEBVs. It must contain `id_col` and all columns listed in
#'   `trait_cols`. Each row should represent one genotype and each trait column
#'   should be numeric or coercible to numeric.
#'
#' @param trait_cols Character vector giving the trait GEBV column names to use
#'   in the index. These names must be present in `gebv_data`. The order of
#'   `trait_cols` determines the order used for `dg`, `W_d`, and any diagonal
#'   weights.
#'
#' @param id_col Character string naming the genotype identifier column used to
#'   align `init_data` and `gebv_data`. Default is `"GenoID"`. This column must
#'   be present in both input tables.
#'
#' @param dg Named numeric desired-gain vector, with names matching
#'   `trait_cols`. This replaces the classical economic-weight vector in the
#'   linear part of the index. For best interpretability, `dg` should be on the
#'   same scale as the transformed GEBV matrix used by the function.
#'
#' @param W_d Optional full symmetric quadratic desired-gain matrix of dimension
#'   `length(trait_cols) x length(trait_cols)`. This matrix controls the
#'   quadratic and interaction contributions among traits. If `NULL`, a diagonal
#'   matrix is created from `quadratic_diag_weights`.
#'
#' @param quadratic_diag_weights Optional named numeric vector of diagonal
#'   weights used only when `W_d = NULL`. Names must match `trait_cols`. Larger
#'   values place stronger emphasis on the quadratic contribution of the
#'   corresponding trait.
#'
#' @param lower_is_better Optional character vector naming traits for which
#'   smaller values are favorable. These traits are internally multiplied by
#'   `-1` so that all traits are aligned in a common favorable direction before
#'   building the index.
#'
#' @param center_traits Logical. If `TRUE`, the GEBV matrix is centered trait by
#'   trait using the column means before the index is computed. This is useful
#'   when GEBVs are not already centered around zero.
#'
#' @param scale_traits Logical. If `TRUE`, each trait column is divided by its
#'   standard deviation after optional centering. This is useful when GEBVs are
#'   on different scales or when you want a standardized trait space.
#'
#' @param impute_missing Logical. If `TRUE`, missing trait values are imputed by
#'   the mean of the corresponding trait column. If `FALSE`, the function stops
#'   with an error when missing values are present.
#'
#' @param return_components Logical. If `TRUE`, the returned ranked table
#'   includes the separate linear and quadratic components in addition to the
#'   final QGSI score. If `FALSE`, only the final QGSI score and its rank are
#'   retained.
#'
#' @param debug Logical. If `TRUE`, the function prints progress and diagnostic
#'   messages during execution. This is useful for tracing data preparation and
#'   matrix construction steps.
#'
#' @return A list containing:
#' \describe{
#'   \item{dg}{The desired-gain vector used in the linear component.}
#'   \item{W_d}{The quadratic desired-gain matrix used in the quadratic component.}
#'   \item{trait_space_used}{A `data.table` containing the transformed GEBV
#'     matrix actually used in the index calculation after optional flipping,
#'     centering, and scaling.}
#'   \item{ranked_geno}{A `data.table` containing `init_data` augmented with
#'     `QGSI_DG`, `Rank_QGSI_DG`, and optionally `LinearDGPart`,
#'     `QuadraticDGPart`, `Rank_LinearDGPart`, and `Rank_QuadraticDGPart`,
#'     sorted by decreasing `QGSI_DG`.}
#'   \item{component_summary}{A `data.table` summarizing the mean, standard
#'     deviation, minimum, and maximum of the computed index components.}
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
#' qdiag <- c(
#'   YLD = 2.0,
#'   MY  = 1.0,
#'   MI  = 1.0,
#'   BL  = 1.5,
#'   NBL = 1.5,
#'   VHB = 1.5
#' )
#'
#' ext <- system.file("extdata", package = "DGQGSI")
#' pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
#' gebv  <- data.table::fread(file.path(ext, "example_gebv.csv"))
#'
#' res <- run_qgsi_desired_gain(
#'   init_data = pheno[, .(GenoID, Family)],
#'   gebv_data = gebv[, c("GenoID", trait_cols), with = FALSE],
#'   trait_cols = trait_cols,
#'   dg = dg,
#'   quadratic_diag_weights = qdiag,
#'   lower_is_better = c("BL", "NBL", "VHB"),
#'   center_traits = FALSE,
#'   scale_traits = FALSE,
#'   impute_missing = TRUE,
#'   return_components = TRUE,
#'   debug = FALSE
#' )
#'
#' head(res$ranked_geno)
#' res$component_summary
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
  .dgqgsi_dbg(debug, "============================================================")
  .dgqgsi_dbg(debug, "Starting run_qgsi_desired_gain()")
  .dgqgsi_dbg(debug, "Traits: %s", paste(trait_cols, collapse = ", "))
  .dgqgsi_dbg(debug, "center_traits = %s | scale_traits = %s", center_traits, scale_traits)
  
  dt_init <- data.table::as.data.table(data.table::copy(init_data))
  dt_gebv <- data.table::as.data.table(data.table::copy(gebv_data))
  
  if (!id_col %in% names(dt_init)) {
    stop(sprintf("id_col '%s' not found in init_data.", id_col), call. = FALSE)
  }
  if (!id_col %in% names(dt_gebv)) {
    stop(sprintf("id_col '%s' not found in gebv_data.", id_col), call. = FALSE)
  }
  if (!all(trait_cols %in% names(dt_gebv))) {
    miss <- setdiff(trait_cols, names(dt_gebv))
    stop("Missing trait columns in gebv_data: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  
  common_ids <- intersect(dt_init[[id_col]], dt_gebv[[id_col]])
  if (length(common_ids) == 0L) {
    stop("No common IDs between init_data and gebv_data.", call. = FALSE)
  }
  
  dt_init <- dt_init[get(id_col) %in% common_ids]
  dt_gebv <- dt_gebv[get(id_col) %in% common_ids]
  data.table::setkeyv(dt_init, id_col)
  data.table::setkeyv(dt_gebv, id_col)
  dt_init <- dt_init[dt_gebv[[id_col]]]
  
  if (!identical(dt_init[[id_col]], dt_gebv[[id_col]])) {
    stop("Alignment failed between init_data and gebv_data.", call. = FALSE)
  }
  
  for (tr in trait_cols) {
    dt_gebv[, (tr) := as.numeric(get(tr))]
  }
  
  if (isTRUE(impute_missing)) {
    for (tr in trait_cols) {
      if (anyNA(dt_gebv[[tr]])) {
        dt_gebv[is.na(get(tr)), (tr) := mean(dt_gebv[[tr]], na.rm = TRUE)]
      }
    }
  } else if (anyNA(as.matrix(dt_gebv[, trait_cols, with = FALSE]))) {
    stop("Missing values detected in gebv_data while impute_missing = FALSE.", call. = FALSE)
  }
  
  if (!is.null(lower_is_better)) {
    if (!all(lower_is_better %in% trait_cols)) {
      bad <- setdiff(lower_is_better, trait_cols)
      stop("These lower_is_better traits are not in trait_cols: ", paste(bad, collapse = ", "), call. = FALSE)
    }
    dt_gebv[, (lower_is_better) := lapply(.SD, function(x) -x), .SDcols = lower_is_better]
  }
  
  gebv_mat <- as.matrix(dt_gebv[, trait_cols, with = FALSE])
  
  if (isTRUE(center_traits)) {
    gebv_mat <- sweep(gebv_mat, 2, colMeans(gebv_mat), FUN = "-")
  }
  
  if (isTRUE(scale_traits)) {
    sds <- apply(gebv_mat, 2, stats::sd)
    bad_sd <- !is.finite(sds) | sds == 0
    sds[bad_sd] <- 1
    gebv_mat <- sweep(gebv_mat, 2, sds, FUN = "/")
  }
  
  colnames(gebv_mat) <- trait_cols
  
  dg <- dg[trait_cols]
  if (anyNA(dg)) {
    stop("dg missing for traits: ", paste(trait_cols[is.na(dg)], collapse = ", "), call. = FALSE)
  }
  dg <- as.numeric(dg)
  names(dg) <- trait_cols
  
  p <- length(trait_cols)
  
  if (is.null(W_d)) {
    if (is.null(quadratic_diag_weights)) {
      quadratic_diag_weights <- rep(1, p)
      names(quadratic_diag_weights) <- trait_cols
    } else {
      quadratic_diag_weights <- quadratic_diag_weights[trait_cols]
      if (anyNA(quadratic_diag_weights)) {
        stop("quadratic_diag_weights missing for traits.", call. = FALSE)
      }
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
      Mean = mean(qgsi_dg),
      SD = stats::sd(qgsi_dg),
      Min = min(qgsi_dg),
      Max = max(qgsi_dg)
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