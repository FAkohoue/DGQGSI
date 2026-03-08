#' Run a full desired-gain index and desired-gain QGSI pipeline
#'
#' @description
#' High-level wrapper that runs one or both of the desired-gain optimization
#' workflow and the desired-gain QGSI workflow, with optional merged comparison output.
#'
#' @param mode One of `"both"`, `"dg"`, or `"qgsi"`.
#' @param init_data A `data.frame` or `data.table` with genotype IDs and metadata.
#' @param trait_cols Character vector of trait names.
#' @param dg Named numeric desired-gain vector.
#' @param id_col Character identifier column. Default is `"GenoID"`.
#' @param cand_data Optional candidate trait data for the DG workflow.
#' @param ref_data Optional reference data for the DG workflow.
#' @param G Optional genetic covariance matrix for the DG workflow.
#' @param select_mode Selection mode for the DG workflow.
#' @param n_select Number selected when `select_mode = "top_n"`.
#' @param trait_min_sd Optional trait thresholds for DG threshold selection.
#' @param fallback_to_top_n Logical; fallback behavior for threshold selection.
#' @param n_iter Iterations for DG optimization.
#' @param n_rep Replicates for DG optimization.
#' @param sd_scale Perturbation scale for DG optimization.
#' @param seed Random seed for DG optimization.
#' @param ridge_P Ridge penalty for DG phenotypic matrix.
#' @param ridge_M Ridge penalty for DG inner matrix.
#' @param dg_scale_traits Logical; if `TRUE`, scale DG input traits.
#' @param gebv_data Optional trait GEBV table for QGSI.
#' @param W_d Optional quadratic desired-gain matrix for QGSI.
#' @param quadratic_diag_weights Optional diagonal quadratic weights for QGSI when `W_d = NULL`.
#' @param qgsi_center_traits Logical; if `TRUE`, center GEBVs before QGSI.
#' @param qgsi_scale_traits Logical; if `TRUE`, scale GEBVs before QGSI.
#' @param qgsi_impute_missing Logical; if `TRUE`, impute missing GEBVs before QGSI.
#' @param W_method Method for automatic `W_d` construction.
#' @param W_base_diag Optional baseline diagonal weights for automatic `W_d`.
#' @param W_lambda_diag Weight for diagonal component in automatic `W_d`.
#' @param W_lambda_outer Weight for outer-product component in automatic `W_d`.
#' @param W_lambda_corr Weight for correlation-weighted component in automatic `W_d`.
#' @param W_corr_power Correlation exponent in automatic `W_d`.
#' @param W_offdiag_shrink Off-diagonal shrinkage in automatic `W_d`.
#' @param W_positive_diagonal_only Logical; if `TRUE`, force non-negative `W_d` diagonal.
#' @param W_normalize Normalization method for automatic `W_d`.
#' @param lower_is_better Optional traits where lower values are preferred.
#' @param merge_outputs Logical; if `TRUE` and `mode = "both"`, merge DG and QGSI outputs.
#' @param compare_sort_by Sorting choice passed to [compare_dg_and_qgsi()].
#' @param debug Logical; if `TRUE`, print debug messages.
#'
#' @return A list containing some or all of `dg_result`, `W_d_result`, `qgsi_result`, and `comparison_result`.
#'
#' @export
run_dgsi_qgsi_pipeline <- function(
    mode = c("both", "dg", "qgsi"),
    init_data,
    trait_cols,
    dg,
    id_col = "GenoID",
    cand_data = NULL,
    ref_data = NULL,
    G = NULL,
    select_mode = c("top_n", "trait_thresholds"),
    n_select = 100,
    trait_min_sd = NULL,
    fallback_to_top_n = TRUE,
    n_iter = 1000,
    n_rep = 20,
    sd_scale = 1.0,
    seed = 42,
    ridge_P = 1e-6,
    ridge_M = 1e-6,
    dg_scale_traits = FALSE,
    gebv_data = NULL,
    W_d = NULL,
    quadratic_diag_weights = NULL,
    qgsi_center_traits = FALSE,
    qgsi_scale_traits = FALSE,
    qgsi_impute_missing = TRUE,
    W_method = c("hybrid", "diag_only", "dg_outer", "corr_weighted"),
    W_base_diag = NULL,
    W_lambda_diag = 1.0,
    W_lambda_outer = 0.5,
    W_lambda_corr = 1.0,
    W_corr_power = 1,
    W_offdiag_shrink = 1.0,
    W_positive_diagonal_only = TRUE,
    W_normalize = c("trace", "max_abs", "none"),
    lower_is_better = NULL,
    merge_outputs = TRUE,
    compare_sort_by = c("DG_rank", "QGSI_rank", "DG", "QGSI", "none"),
    debug = TRUE
) {
  mode <- match.arg(mode)
  select_mode <- match.arg(select_mode)
  W_method <- match.arg(W_method)
  W_normalize <- match.arg(W_normalize)
  compare_sort_by <- match.arg(compare_sort_by)

  out <- list()

  if (mode %in% c("dg", "both")) {
    if (is.null(cand_data)) stop("cand_data must be provided when mode includes 'dg'.", call. = FALSE)
    out$dg_result <- run_desired_gain_index_Joukhadar2024(
      init_data = init_data,
      cand_data = cand_data,
      trait_cols = trait_cols,
      ref_data = ref_data,
      id_col = id_col,
      scale_traits = dg_scale_traits,
      lower_is_better = lower_is_better,
      G = G,
      dg = dg,
      select_mode = select_mode,
      n_select = n_select,
      trait_min_sd = trait_min_sd,
      fallback_to_top_n = fallback_to_top_n,
      n_iter = n_iter,
      n_rep = n_rep,
      sd_scale = sd_scale,
      seed = seed,
      ridge_P = ridge_P,
      ridge_M = ridge_M,
      debug = debug,
      return_all_reps = TRUE
    )
  }

  if (mode %in% c("qgsi", "both")) {
    if (is.null(gebv_data)) stop("gebv_data must be provided when mode includes 'qgsi'.", call. = FALSE)
    if (is.null(W_d)) {
      out$W_d_result <- construct_Wd_from_dg(
        gebv_data = gebv_data,
        trait_cols = trait_cols,
        dg = dg,
        lower_is_better = lower_is_better,
        center_traits = qgsi_center_traits,
        scale_traits = qgsi_scale_traits,
        impute_missing = qgsi_impute_missing,
        method = W_method,
        base_diag = W_base_diag,
        lambda_diag = W_lambda_diag,
        lambda_outer = W_lambda_outer,
        lambda_corr = W_lambda_corr,
        corr_power = W_corr_power,
        offdiag_shrink = W_offdiag_shrink,
        positive_diagonal_only = W_positive_diagonal_only,
        normalize = W_normalize,
        debug = debug
      )
      W_d <- out$W_d_result$W_d
    }

    out$qgsi_result <- run_qgsi_desired_gain(
      init_data = init_data,
      gebv_data = gebv_data,
      trait_cols = trait_cols,
      id_col = id_col,
      dg = dg,
      W_d = W_d,
      quadratic_diag_weights = quadratic_diag_weights,
      lower_is_better = lower_is_better,
      center_traits = qgsi_center_traits,
      scale_traits = qgsi_scale_traits,
      impute_missing = qgsi_impute_missing,
      return_components = TRUE,
      debug = debug
    )
  }

  if (mode == "both" && isTRUE(merge_outputs)) {
    out$comparison_result <- compare_dg_and_qgsi(
      dg_result = out$dg_result,
      qgsi_result = out$qgsi_result,
      id_col = id_col,
      sort_by = compare_sort_by,
      debug = debug
    )
  }

  out
}
