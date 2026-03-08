#' Merge desired-gain index and desired-gain QGSI outputs into one comparison table
#'
#' @description
#' Merges the outputs of [run_desired_gain_index_Joukhadar2024()] and
#' [run_qgsi_desired_gain()] into one genotype-level comparison table by ID.
#'
#' @param dg_result Result list returned by `run_desired_gain_index_Joukhadar2024()`.
#' @param qgsi_result Result list returned by `run_qgsi_desired_gain()`.
#' @param id_col Character identifier column. Default is `"GenoID"`.
#' @param dg_cols Columns to extract from `dg_result$ranked_geno`.
#' @param qgsi_cols Columns to extract from `qgsi_result$ranked_geno`.
#' @param include_metadata Logical; retained for API consistency.
#' @param compute_rank_differences Logical; if `TRUE`, compute rank differences.
#' @param add_correlation_summary Logical; if `TRUE`, compute score/rank correlations.
#' @param sort_by One of `"DG_rank"`, `"QGSI_rank"`, `"DG"`, `"QGSI"`, or `"none"`.
#' @param debug Logical; if `TRUE`, print debug messages.
#'
#' @return A list with `comparison_table` and optional `correlation_summary`.
#'
#' @import data.table
#' @export
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
