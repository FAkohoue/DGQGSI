\name{run_desired_gain_index_Joukhadar2024}
\alias{run_desired_gain_index_Joukhadar2024}
\title{Run a desired-gain selection index}
\usage{
run_desired_gain_index_Joukhadar2024(init_data, cand_data, trait_cols, ref_data = NULL,
  id_col = "GenoID", scale_traits = FALSE, lower_is_better = NULL,
  G = NULL, dg, select_mode = c("top_n", "trait_thresholds"),
  n_select = 100, trait_min_sd = NULL, fallback_to_top_n = TRUE,
  n_iter = 1000, n_rep = 20, sd_scale = 1, seed = 42,
  ridge_P = 1e-06, ridge_M = 1e-06, debug = TRUE, return_all_reps = TRUE)
}
\arguments{
\item{init_data}{Metadata table.}
\item{cand_data}{Candidate trait table.}
\item{trait_cols}{Trait columns.}
\item{ref_data}{Reference table.}
\item{id_col}{Identifier column.}
\item{scale_traits}{Scale traits or not.}
\item{lower_is_better}{Traits to flip.}
\item{G}{Genetic covariance matrix.}
\item{dg}{Desired-gain vector.}
\item{select_mode}{Selection mode.}
\item{n_select}{Number selected in top-n mode.}
\item{trait_min_sd}{Thresholds in threshold mode.}
\item{fallback_to_top_n}{Fallback behavior.}
\item{n_iter}{Iterations.}
\item{n_rep}{Replicates.}
\item{sd_scale}{Perturbation scale.}
\item{seed}{Random seed.}
\item{ridge_P}{Ridge penalty for P.}
\item{ridge_M}{Ridge penalty for M.}
\item{debug}{Print debug messages.}
\item{return_all_reps}{Return all replicates.}
}
\value{
A list with optimized coefficients, realized gains, ranked output, and optional replicate details.
}
\description{
Optimizes desired-gain index coefficients and returns a ranked genotype table plus selected subsets.
}
\keyword{methods}
