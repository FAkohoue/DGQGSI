\name{construct_Wd_from_dg}
\alias{construct_Wd_from_dg}
\title{Construct a quadratic desired-gain matrix for QGSI}
\usage{
construct_Wd_from_dg(gebv_data, trait_cols, dg, lower_is_better = NULL,
  center_traits = TRUE, scale_traits = TRUE, impute_missing = TRUE,
  method = c("hybrid", "diag_only", "dg_outer", "corr_weighted"),
  base_diag = NULL, lambda_diag = 1, lambda_outer = 0.5,
  lambda_corr = 1, corr_power = 1, offdiag_shrink = 1,
  positive_diagonal_only = TRUE, normalize = c("trace", "max_abs", "none"),
  debug = TRUE)
}
\arguments{
\item{gebv_data}{Trait GEBV table.}
\item{trait_cols}{Trait columns.}
\item{dg}{Desired-gain vector.}
\item{lower_is_better}{Traits to flip.}
\item{center_traits}{Center traits.}
\item{scale_traits}{Scale traits.}
\item{impute_missing}{Impute missing values.}
\item{method}{Construction method.}
\item{base_diag}{Optional base diagonal.}
\item{lambda_diag}{Diagonal weight.}
\item{lambda_outer}{Outer-product weight.}
\item{lambda_corr}{Correlation-weighted component weight.}
\item{corr_power}{Correlation exponent.}
\item{offdiag_shrink}{Off-diagonal shrinkage.}
\item{positive_diagonal_only}{Force non-negative diagonal.}
\item{normalize}{Normalization method.}
\item{debug}{Print debug messages.}
}
\value{
A list containing W_d and intermediate components.
}
\description{
Builds a symmetric quadratic desired-gain matrix for the desired-gain QGSI formulation.
}
\keyword{methods}
