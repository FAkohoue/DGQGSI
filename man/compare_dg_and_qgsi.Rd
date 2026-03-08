\name{compare_dg_and_qgsi}
\alias{compare_dg_and_qgsi}
\title{Compare DG and QGSI outputs}
\usage{
compare_dg_and_qgsi(dg_result, qgsi_result, id_col = "GenoID",
  dg_cols = c(id_col, "SelectionIndex", "Selected"),
  qgsi_cols = c(id_col, "LinearDGPart", "QuadraticDGPart", "QGSI_DG", "Rank_QGSI_DG"),
  include_metadata = TRUE, compute_rank_differences = TRUE,
  add_correlation_summary = TRUE, sort_by = c("DG_rank", "QGSI_rank", "DG", "QGSI", "none"),
  debug = TRUE)
}
\arguments{
\item{dg_result}{Output of run_desired_gain_index_Joukhadar2024().}
\item{qgsi_result}{Output of run_qgsi_desired_gain().}
\item{id_col}{Identifier column.}
\item{dg_cols}{Columns to extract from DG result.}
\item{qgsi_cols}{Columns to extract from QGSI result.}
\item{include_metadata}{Retained for API consistency.}
\item{compute_rank_differences}{Compute rank differences.}
\item{add_correlation_summary}{Compute score/rank correlations.}
\item{sort_by}{Final sort order.}
\item{debug}{Print debug messages.}
}
\value{
A list with comparison_table and optional correlation_summary.
}
\description{
Merges DG and QGSI outputs into a single comparison table by genotype ID.
}
\keyword{methods}
