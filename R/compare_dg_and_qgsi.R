#' Merge desired-gain index and desired-gain QGSI outputs into one comparison table
#'
#' @description
#' Merges the outputs of [run_desired_gain_index_Joukhadar2024()] and
#' [run_qgsi_desired_gain()] into one genotype-level comparison table by ID.
#'
#' This helper function is intended for situations where both the
#' desired-gain selection index workflow and the desired-gain QGSI workflow
#' have been run on the same or overlapping sets of genotypes, and the user
#' wants to compare:
#' \itemize{
#'   \item genotype scores from both methods,
#'   \item genotype ranks from both methods,
#'   \item selected versus non-selected status from the DGSI workflow,
#'   \item agreement or disagreement between ranking systems.
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item extracts selected columns from `dg_result$ranked_geno`;
#'   \item extracts selected columns from `qgsi_result$ranked_geno`;
#'   \item renames overlapping DG columns for clarity;
#'   \item merges both tables by `id_col`;
#'   \item optionally computes rank-difference diagnostics;
#'   \item optionally computes correlation summaries between DGSI and QGSI scores
#'         and/or ranks;
#'   \item optionally sorts the final merged table according to the chosen
#'         ranking criterion.
#' }
#'
#' This is especially useful for downstream interpretation, because it allows
#' side-by-side comparison of:
#' \itemize{
#'   \item `DG_SelectionIndex`
#'   \item `QGSI_DG`
#'   \item their respective ranks
#'   \item and the magnitude of rank disagreement
#' }
#'
#' If some requested columns are not found in either result object, they are
#' skipped with a debug message rather than causing immediate failure, as long
#' as the identifier column is available.
#'
#' @param dg_result A result list returned by
#'   [run_desired_gain_index_Joukhadar2024()]. The function expects at least a
#'   `ranked_geno` component containing the DGSI output table.
#'
#' @param qgsi_result A result list returned by [run_qgsi_desired_gain()].
#'   The function expects at least a `ranked_geno` component containing the QGSI
#'   output table.
#'
#' @param id_col Character string naming the genotype identifier column used to
#'   merge DGSI and QGSI outputs. Default is `"GenoID"`. This column must be
#'   present in both `dg_result$ranked_geno` and `qgsi_result$ranked_geno`.
#'
#' @param dg_cols Character vector giving the columns to extract from
#'   `dg_result$ranked_geno`. By default, this includes the genotype identifier,
#'   the DG selection index, and the selected/non-selected flag. Additional
#'   metadata columns such as family or population can also be requested here.
#'
#' @param qgsi_cols Character vector giving the columns to extract from
#'   `qgsi_result$ranked_geno`. By default, this includes the genotype
#'   identifier, the QGSI linear and quadratic components, the final QGSI score,
#'   and the QGSI rank.
#'
#' @param include_metadata Logical retained for API consistency and future use.
#'   At present, metadata are included whenever they are explicitly requested in
#'   `dg_cols` or `qgsi_cols`. Default is `TRUE`.
#'
#' @param compute_rank_differences Logical. If `TRUE`, the function computes
#'   diagnostic columns comparing DGSI and QGSI ranks, including:
#'   \itemize{
#'     \item `RankDiff_DG_minus_QGSI`
#'     \item `AbsRankDiff_DG_vs_QGSI`
#'     \item `SignSame_DG_QGSI` when both score columns are available
#'   }
#'
#' @param add_correlation_summary Logical. If `TRUE`, the function computes a
#'   summary table containing Pearson and Spearman correlations between:
#'   \itemize{
#'     \item DG selection index and QGSI score
#'     \item DG rank and QGSI rank
#'   }
#'   whenever the required columns are available.
#'
#' @param sort_by Character string specifying how the merged output table should
#'   be ordered. One of:
#'   \itemize{
#'     \item `"DG_rank"`: sort by DG rank ascending
#'     \item `"QGSI_rank"`: sort by QGSI rank ascending
#'     \item `"DG"`: sort by DG score descending
#'     \item `"QGSI"`: sort by QGSI score descending
#'     \item `"none"`: do not apply additional sorting after merging
#'   }
#'
#' @param debug Logical. If `TRUE`, the function prints progress and diagnostic
#'   messages, including missing requested columns, computed correlation
#'   summaries, and the final sorting rule applied.
#'
#' @return A list containing:
#' \describe{
#'   \item{comparison_table}{A `data.table` containing the merged genotype-level
#'     outputs from DGSI and QGSI, optionally augmented with rank-difference
#'     diagnostics.}
#'   \item{correlation_summary}{A `data.table` containing Pearson and Spearman
#'     correlations between DGSI and QGSI scores and/or ranks when
#'     `add_correlation_summary = TRUE`; otherwise `NULL`.}
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
#' trait_min_sd <- c(
#'   YLD = 0.2,
#'   MY  = 0.1,
#'   MI  = 0.1,
#'   BL  = 0.1,
#'   NBL = 0.1,
#'   VHB = 0.1
#' )
#'
#' ext <- system.file("extdata", package = "DGQGSI")
#' pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
#' gebv  <- data.table::fread(file.path(ext, "example_gebv.csv"))
#'
#' dg_res <- run_desired_gain_index_Joukhadar2024(
#'   init_data = pheno[, .(GenoID, Family)],
#'   cand_data = pheno[, c("GenoID", trait_cols), with = FALSE],
#'   ref_data = pheno[, c("GenoID", trait_cols), with = FALSE],
#'   trait_cols = trait_cols,
#'   dg = dg,
#'   lower_is_better = c("BL", "NBL", "VHB"),
#'   select_mode = "trait_thresholds",
#'   trait_min_sd = trait_min_sd,
#'   scale_traits = FALSE,
#'   debug = FALSE
#' )
#'
#' qgsi_res <- run_qgsi_desired_gain(
#'   init_data = pheno[, .(GenoID, Family)],
#'   gebv_data = gebv[, c("GenoID", trait_cols), with = FALSE],
#'   trait_cols = trait_cols,
#'   dg = dg,
#'   lower_is_better = c("BL", "NBL", "VHB"),
#'   center_traits = FALSE,
#'   scale_traits = FALSE,
#'   debug = FALSE
#' )
#'
#' merged <- compare_dg_and_qgsi(
#'   dg_result = dg_res,
#'   qgsi_result = qgsi_res,
#'   id_col = "GenoID",
#'   compute_rank_differences = TRUE,
#'   add_correlation_summary = TRUE,
#'   sort_by = "DG_rank",
#'   debug = FALSE
#' )
#'
#' head(merged$comparison_table)
#' merged$correlation_summary
#'
#' @seealso
#' [run_desired_gain_index_Joukhadar2024()],
#' [run_qgsi_desired_gain()],
#' [run_dgsi_qgsi_pipeline()]
#'
#' @import data.table
#' @export
#' 
compare_dg_and_qgsi <- function(
    dg_result,
    qgsi_result,
    id_col = "GenoID",
    dg_cols = c(id_col, "SelectionIndex", "Selected"),
    qgsi_cols = c(id_col, "LinearDGPart", "QuadraticDGPart", "QGSI_DG", "Rank_QGSI_DG"),
    include_metadata = TRUE,
    compute_rank_differences = TRUE,
    add_correlation_summary = TRUE,
    sort_by = c("DG_rank", "QGSI_rank", "DG", "QGSI", "none"),
    debug = TRUE
) {
  sort_by <- match.arg(sort_by)

  if (is.null(dg_result$ranked_geno)) stop("dg_result does not contain 'ranked_geno'.", call. = FALSE)
  if (is.null(qgsi_result$ranked_geno)) stop("qgsi_result does not contain 'ranked_geno'.", call. = FALSE)

  dt_dg <- data.table::as.data.table(data.table::copy(dg_result$ranked_geno))
  dt_qg <- data.table::as.data.table(data.table::copy(qgsi_result$ranked_geno))

  if (!id_col %in% names(dt_dg)) stop(sprintf("id_col '%s' not found in dg_result$ranked_geno.", id_col), call. = FALSE)
  if (!id_col %in% names(dt_qg)) stop(sprintf("id_col '%s' not found in qgsi_result$ranked_geno.", id_col), call. = FALSE)

  if (!"Rank_SelectionIndex" %in% names(dt_dg) && "SelectionIndex" %in% names(dt_dg)) {
    dt_dg[, Rank_SelectionIndex := data.table::frank(-SelectionIndex, ties.method = "average")]
  }

  dg_cols_present <- intersect(dg_cols, names(dt_dg))
  qgsi_cols_present <- intersect(qgsi_cols, names(dt_qg))
  if (!id_col %in% dg_cols_present) dg_cols_present <- unique(c(id_col, dg_cols_present))
  if (!id_col %in% qgsi_cols_present) qgsi_cols_present <- unique(c(id_col, qgsi_cols_present))

  dg_sub <- dt_dg[, ..dg_cols_present]
  qgsi_sub <- dt_qg[, ..qgsi_cols_present]

  if ("SelectionIndex" %in% names(dg_sub)) data.table::setnames(dg_sub, "SelectionIndex", "DG_SelectionIndex")
  if ("Rank_SelectionIndex" %in% names(dg_sub)) data.table::setnames(dg_sub, "Rank_SelectionIndex", "DG_Rank_SelectionIndex")
  if ("Selected" %in% names(dg_sub)) data.table::setnames(dg_sub, "Selected", "DG_Selected")

  merged_dt <- merge(dg_sub, qgsi_sub, by = id_col, all = TRUE, sort = FALSE)

  if (isTRUE(compute_rank_differences)) {
    if ("DG_Rank_SelectionIndex" %in% names(merged_dt) && "Rank_QGSI_DG" %in% names(merged_dt)) {
      merged_dt[, RankDiff_DG_minus_QGSI := DG_Rank_SelectionIndex - Rank_QGSI_DG]
      merged_dt[, AbsRankDiff_DG_vs_QGSI := abs(RankDiff_DG_minus_QGSI)]
    }
    if ("DG_SelectionIndex" %in% names(merged_dt) && "QGSI_DG" %in% names(merged_dt)) {
      merged_dt[, SignSame_DG_QGSI := sign(DG_SelectionIndex) == sign(QGSI_DG)]
    }
  }

  correlation_summary <- NULL
  if (isTRUE(add_correlation_summary)) {
    corr_rows <- list()
    if ("DG_SelectionIndex" %in% names(merged_dt) && "QGSI_DG" %in% names(merged_dt)) {
      x <- merged_dt$DG_SelectionIndex; y <- merged_dt$QGSI_DG
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) >= 3) {
        corr_rows[[length(corr_rows)+1L]] <- data.table::data.table(
          Comparison = "DG_SelectionIndex vs QGSI_DG",
          Pearson = stats::cor(x[ok], y[ok], method = "pearson"),
          Spearman = stats::cor(x[ok], y[ok], method = "spearman"),
          N = sum(ok)
        )
      }
    }
    if ("DG_Rank_SelectionIndex" %in% names(merged_dt) && "Rank_QGSI_DG" %in% names(merged_dt)) {
      x <- merged_dt$DG_Rank_SelectionIndex; y <- merged_dt$Rank_QGSI_DG
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) >= 3) {
        corr_rows[[length(corr_rows)+1L]] <- data.table::data.table(
          Comparison = "DG_Rank_SelectionIndex vs Rank_QGSI_DG",
          Pearson = stats::cor(x[ok], y[ok], method = "pearson"),
          Spearman = stats::cor(x[ok], y[ok], method = "spearman"),
          N = sum(ok)
        )
      }
    }
    if (length(corr_rows) > 0L) correlation_summary <- data.table::rbindlist(corr_rows, fill = TRUE)
  }

  if (sort_by == "DG" && "DG_SelectionIndex" %in% names(merged_dt)) data.table::setorder(merged_dt, -DG_SelectionIndex)
  if (sort_by == "QGSI" && "QGSI_DG" %in% names(merged_dt)) data.table::setorder(merged_dt, -QGSI_DG)
  if (sort_by == "DG_rank" && "DG_Rank_SelectionIndex" %in% names(merged_dt)) data.table::setorder(merged_dt, DG_Rank_SelectionIndex)
  if (sort_by == "QGSI_rank" && "Rank_QGSI_DG" %in% names(merged_dt)) data.table::setorder(merged_dt, Rank_QGSI_DG)

  list(comparison_table = merged_dt, correlation_summary = correlation_summary)
}
