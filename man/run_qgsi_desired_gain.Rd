\name{run_qgsi_desired_gain}
\alias{run_qgsi_desired_gain}
\title{Run desired-gain QGSI using trait GEBVs}
\usage{
run_qgsi_desired_gain(init_data, gebv_data, trait_cols, id_col = "GenoID",
  dg, W_d = NULL, quadratic_diag_weights = NULL, lower_is_better = NULL,
  center_traits = TRUE, scale_traits = TRUE, impute_missing = TRUE,
  return_components = TRUE, debug = TRUE)
}
\arguments{
\item{init_data}{Metadata table.}
\item{gebv_data}{Trait GEBV table.}
\item{trait_cols}{Trait columns.}
\item{id_col}{Identifier column.}
\item{dg}{Desired-gain vector.}
\item{W_d}{Quadratic matrix.}
\item{quadratic_diag_weights}{Diagonal weights if W_d is NULL.}
\item{lower_is_better}{Traits to flip.}
\item{center_traits}{Center GEBVs.}
\item{scale_traits}{Scale GEBVs.}
\item{impute_missing}{Impute missing GEBVs.}
\item{return_components}{Return linear and quadratic parts.}
\item{debug}{Print debug messages.}
}
\value{
A list containing the transformed trait space, ranked genotypes, and component summaries.
}
\description{
Computes a desired-gain adaptation of the quadratic genomic selection index.
}
\keyword{methods}
