#' Run a full desired-gain index and desired-gain QGSI pipeline
#'
#' @description
#' High-level wrapper that runs one or both of the desired-gain selection index
#' optimization workflow and the desired-gain quadratic genomic selection index
#' (QGSI) workflow, with optional merged comparison output.
#'
#' This function orchestrates the full DGQGSI pipeline, allowing users to:
#' \itemize{
#'   \item compute a desired-gain selection index following
#'   Joukhadar et al. (2024),
#'   \item compute a desired-gain quadratic genomic selection index (QGSI)
#'   following Cerón-Rojas et al. (2026),
#'   \item automatically construct the quadratic desired-gain matrix `W_d`
#'   when needed,
#'   \item and generate a merged comparison table between the two methods.
#' }
#'
#' The function acts as a convenient wrapper that coordinates multiple package
#' functions, ensuring consistent input alignment and reducing repetitive code
#' when running complete selection-index analyses.
#'
#' @details
#' The pipeline can operate in three modes controlled by `mode`:
#'
#' \describe{
#' \item{`"dg"`}{
#' Runs only the desired-gain selection index workflow using
#' [run_desired_gain_index_Joukhadar2024()].
#' }
#'
#' \item{`"qgsi"`}{
#' Runs only the desired-gain QGSI workflow using
#' [run_qgsi_desired_gain()]. If `W_d` is not provided, the matrix is
#' automatically constructed using [construct_Wd_from_dg()].
#' }
#'
#' \item{`"both"`}{
#' Runs both workflows sequentially. If `merge_outputs = TRUE`,
#' the results are merged using [compare_dg_and_qgsi()] to produce a
#' genotype-level comparison table.
#' }
#' }
#'
#' The typical pipeline sequence when `mode = "both"` is:
#'
#' \enumerate{
#'   \item run the desired-gain index optimization;
#'   \item optionally construct a quadratic desired-gain matrix `W_d`;
#'   \item compute the QGSI scores using GEBVs;
#'   \item optionally merge both outputs into a comparison table.
#' }
#'
#' This wrapper is particularly useful in genomic selection pipelines
#' where breeders want to compare:
#' \itemize{
#'   \item classical linear desired-gain indices,
#'   \item nonlinear quadratic genomic selection indices,
#'   \item and their resulting genotype rankings.
#' }
#'
#' @param mode Character string specifying which workflow to run.
#'   One of `"both"`, `"dg"`, or `"qgsi"`.
#'
#' @param init_data A `data.frame` or `data.table` containing genotype
#'   identifiers and metadata (for example `GenoID`, family, population).
#'   This table is preserved and appended with index results.
#'
#' @param trait_cols Character vector of trait names used in both DG and
#'   QGSI calculations. These columns must exist in the candidate trait
#'   data (`cand_data`) and/or the GEBV data (`gebv_data`).
#'
#' @param dg Named numeric desired-gain vector whose names correspond to
#'   `trait_cols`. This vector defines the breeding objectives and is used
#'   in both DGSI and QGSI calculations.
#'
#' @param id_col Character genotype identifier column used to align all
#'   input tables. Default is `"GenoID"`.
#'
#' @param cand_data Optional `data.frame` or `data.table` containing
#'   candidate trait values for the desired-gain index workflow.
#'
#' @param ref_data Optional reference dataset used for scaling and
#'   covariance estimation in the desired-gain index workflow.
#'
#' @param G Optional genetic covariance matrix used in the desired-gain
#'   index optimization. If `NULL`, it is estimated from the reference data.
#'
#' @param select_mode Selection strategy used in the DG workflow.
#'   One of `"top_n"` or `"trait_thresholds"`.
#'
#' @param n_select Integer specifying the number of selected genotypes
#'   when `select_mode = "top_n"`.
#'
#' @param trait_min_sd Named numeric vector specifying minimum selection
#'   thresholds when `select_mode = "trait_thresholds"`. Values are
#'   expressed in standard deviation units.
#'
#' @param fallback_to_top_n Logical. If `TRUE`, the function falls back to
#'   top-N selection when threshold selection returns zero genotypes.
#'
#' @param n_iter Integer number of optimization iterations used in the
#'   desired-gain index search algorithm.
#'
#' @param n_rep Integer number of optimization replicates.
#'
#' @param sd_scale Numeric perturbation scale applied to candidate
#'   desired-gain vectors during optimization.
#'
#' @param seed Integer random seed controlling reproducibility of the
#'   optimization procedure.
#'
#' @param ridge_P Numeric ridge penalty added to the phenotypic correlation
#'   matrix used in the desired-gain index calculation.
#'
#' @param ridge_M Numeric ridge penalty added to the inner matrix used
#'   during coefficient computation.
#'
#' @param dg_scale_traits Logical. If `TRUE`, traits are centered and
#'   scaled before computing the desired-gain index.
#'
#' @param gebv_data Optional `data.frame` or `data.table` containing
#'   genotype-level GEBVs used for QGSI computation.
#'
#' @param W_d Optional quadratic desired-gain matrix. If `NULL`, the
#'   matrix is automatically generated using
#'   [construct_Wd_from_dg()].
#'
#' @param quadratic_diag_weights Optional named numeric vector used to
#'   construct a diagonal `W_d` matrix when no full matrix is provided.
#'
#' @param qgsi_center_traits Logical. If `TRUE`, the GEBV matrix is
#'   centered before computing QGSI.
#'
#' @param qgsi_scale_traits Logical. If `TRUE`, the GEBV matrix is
#'   scaled before computing QGSI.
#'
#' @param qgsi_impute_missing Logical. If `TRUE`, missing GEBV values
#'   are imputed by trait means before computing QGSI.
#'
#' @param W_method Method used to automatically construct `W_d`.
#'   Passed to [construct_Wd_from_dg()]. Options include
#'   `"hybrid"`, `"diag_only"`, `"dg_outer"`, and `"corr_weighted"`.
#'
#' @param W_base_diag Optional baseline diagonal weights for automatic
#'   `W_d` construction.
#'
#' @param W_lambda_diag Weight applied to the diagonal component when
#'   building `W_d`.
#'
#' @param W_lambda_outer Weight applied to the desired-gain outer-product
#'   component of `W_d`.
#'
#' @param W_lambda_corr Weight applied to the correlation-weighted
#'   component of `W_d`.
#'
#' @param W_corr_power Exponent applied to trait correlations when
#'   constructing the correlation-weighted component.
#'
#' @param W_offdiag_shrink Numeric shrinkage factor applied to
#'   off-diagonal elements of `W_d`.
#'
#' @param W_positive_diagonal_only Logical. If `TRUE`, negative diagonal
#'   entries in `W_d` are set to zero.
#'
#' @param W_normalize Normalization strategy applied to `W_d`. One of
#'   `"trace"`, `"max_abs"`, or `"none"`.
#'
#' @param lower_is_better Optional character vector listing traits for
#'   which smaller values are favorable. These traits are internally
#'   flipped so that all objectives align in the same direction.
#'
#' @param merge_outputs Logical. If `TRUE` and `mode = "both"`, the
#'   DG and QGSI outputs are merged using
#'   [compare_dg_and_qgsi()].
#'
#' @param compare_sort_by Sorting rule passed to
#'   [compare_dg_and_qgsi()]. Determines how the merged comparison table
#'   is ordered.
#'
#' @param debug Logical. If `TRUE`, the function prints progress and
#'   diagnostic messages during execution.
#'
#' @return A list that may contain:
#' \describe{
#'   \item{dg_result}{Output of
#'     [run_desired_gain_index_Joukhadar2024()].}
#'
#'   \item{W_d_result}{Output of
#'     [construct_Wd_from_dg()] when the quadratic matrix is
#'     constructed automatically.}
#'
#'   \item{qgsi_result}{Output of
#'     [run_qgsi_desired_gain()].}
#'
#'   \item{comparison_result}{Output of
#'     [compare_dg_and_qgsi()] when both workflows are run
#'     and `merge_outputs = TRUE`.}
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
#' pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
#' gebv  <- data.table::fread(file.path(ext, "example_gebv.csv"))
#'
#' res <- run_dgsi_qgsi_pipeline(
#'   mode = "both",
#'   init_data = pheno[, .(GenoID, Family)],
#'   cand_data = pheno[, c("GenoID", trait_cols), with = FALSE],
#'   gebv_data = gebv[, c("GenoID", trait_cols), with = FALSE],
#'   trait_cols = trait_cols,
#'   dg = dg,
#'   lower_is_better = c("BL", "NBL", "VHB"),
#'   merge_outputs = TRUE,
#'   debug = FALSE
#' )
#'
#' head(res$comparison_result$comparison_table)
#'
#' @references
#' Joukhadar R et al. (2024).
#' Optimising desired gain indices to maximise selection response.
#' Frontiers in Plant Science 15:1337388.
#'
#' Cerón-Rojas JJ et al. (2026).
#' Nonlinear genomic selection index accelerates multi-trait crop improvement.
#' Nature Communications 17:1991.
#'
#' @seealso
#' [run_desired_gain_index_Joukhadar2024()],
#' [run_qgsi_desired_gain()],
#' [construct_Wd_from_dg()],
#' [compare_dg_and_qgsi()]
#'
#' @export
#' 
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
