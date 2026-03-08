\name{run_dgsi_qgsi_pipeline}
\alias{run_dgsi_qgsi_pipeline}
\title{Run a full DG and QGSI pipeline}
\usage{
run_dgsi_qgsi_pipeline(mode = c("both", "dg", "qgsi"), init_data,
  trait_cols, dg, id_col = "GenoID", cand_data = NULL, ref_data = NULL,
  G = NULL, select_mode = c("top_n", "trait_thresholds"),
  n_select = 100, trait_min_sd = NULL, fallback_to_top_n = TRUE,
  n_iter = 1000, n_rep = 20, sd_scale = 1, seed = 42,
  ridge_P = 1e-06, ridge_M = 1e-06, dg_scale_traits = FALSE,
  gebv_data = NULL, W_d = NULL, quadratic_diag_weights = NULL,
  qgsi_center_traits = FALSE, qgsi_scale_traits = FALSE,
  qgsi_impute_missing = TRUE, W_method = c("hybrid", "diag_only", "dg_outer", "corr_weighted"),
  W_base_diag = NULL, W_lambda_diag = 1, W_lambda_outer = 0.5,
  W_lambda_corr = 1, W_corr_power = 1, W_offdiag_shrink = 1,
  W_positive_diagonal_only = TRUE, W_normalize = c("trace", "max_abs", "none"),
  lower_is_better = NULL, merge_outputs = TRUE,
  compare_sort_by = c("DG_rank", "QGSI_rank", "DG", "QGSI", "none"),
  debug = TRUE)
}
\arguments{
\item{mode}{Run DG only, QGSI only, or both.}
\item{init_data}{Metadata table.}
\item{trait_cols}{Trait columns.}
\item{dg}{Desired-gain vector.}
\item{id_col}{Identifier column.}
\item{cand_data}{Candidate trait table for DG.}
\item{ref_data}{Reference table for DG.}
\item{G}{Genetic covariance matrix for DG.}
\item{select_mode}{Selection mode for DG.}
\item{n_select}{Number selected in top-n mode.}
\item{trait_min_sd}{Thresholds for DG.}
\item{fallback_to_top_n}{Fallback behavior.}
\item{n_iter}{DG iterations.}
\item{n_rep}{DG replicates.}
\item{sd_scale}{DG perturbation scale.}
\item{seed}{Random seed.}
\item{ridge_P}{DG ridge for P.}
\item{ridge_M}{DG ridge for M.}
\item{dg_scale_traits}{Scale DG inputs.}
\item{gebv_data}{Trait GEBV table for QGSI.}
\item{W_d}{Optional quadratic matrix.}
\item{quadratic_diag_weights}{Optional diagonal weights.}
\item{qgsi_center_traits}{Center GEBVs.}
\item{qgsi_scale_traits}{Scale GEBVs.}
\item{qgsi_impute_missing}{Impute missing GEBVs.}
\item{W_method}{Automatic W_d method.}
\item{W_base_diag}{Optional base diagonal.}
\item{W_lambda_diag}{Diagonal weight.}
\item{W_lambda_outer}{Outer-product weight.}
\item{W_lambda_corr}{Correlation-weighted component weight.}
\item{W_corr_power}{Correlation exponent.}
\item{W_offdiag_shrink}{Off-diagonal shrinkage.}
\item{W_positive_diagonal_only}{Force non-negative diagonal.}
\item{W_normalize}{Normalization method.}
\item{lower_is_better}{Traits to flip.}
\item{merge_outputs}{Merge outputs when mode is both.}
\item{compare_sort_by}{Sort order for merged output.}
\item{debug}{Print debug messages.}
}
\value{
A list containing dg_result, W_d_result, qgsi_result, and optionally comparison_result.
}
\description{
High-level package wrapper for one-stop execution of the DG and desired-gain QGSI workflows.
}
\keyword{methods}
