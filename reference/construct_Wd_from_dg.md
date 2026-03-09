# Construct a quadratic desired-gain matrix for QGSI

Builds a symmetric quadratic desired-gain matrix for the desired-gain
QGSI formulation.

## Usage

``` r
construct_Wd_from_dg(gebv_data, trait_cols, dg, lower_is_better = NULL,
  center_traits = TRUE, scale_traits = TRUE, impute_missing = TRUE,
  method = c("hybrid", "diag_only", "dg_outer", "corr_weighted"),
  base_diag = NULL, lambda_diag = 1, lambda_outer = 0.5,
  lambda_corr = 1, corr_power = 1, offdiag_shrink = 1,
  positive_diagonal_only = TRUE, normalize = c("trace", "max_abs", "none"),
  debug = TRUE)
```

## Arguments

- gebv_data:

  Trait GEBV table.

- trait_cols:

  Trait columns.

- dg:

  Desired-gain vector.

- lower_is_better:

  Traits to flip.

- center_traits:

  Center traits.

- scale_traits:

  Scale traits.

- impute_missing:

  Impute missing values.

- method:

  Construction method.

- base_diag:

  Optional base diagonal.

- lambda_diag:

  Diagonal weight.

- lambda_outer:

  Outer-product weight.

- lambda_corr:

  Correlation-weighted component weight.

- corr_power:

  Correlation exponent.

- offdiag_shrink:

  Off-diagonal shrinkage.

- positive_diagonal_only:

  Force non-negative diagonal.

- normalize:

  Normalization method.

- debug:

  Print debug messages.

## Value

A list containing W_d and intermediate components.
