# Construct a quadratic desired-gain matrix `W_d` automatically for QGSI

Builds a symmetric quadratic desired-gain matrix for the desired-gain
QGSI formulation: \$\$QGSI\\DG_i = d' \hat{\gamma}\_i + \hat{\gamma}\_i'
W_d \hat{\gamma}\_i\$\$

This helper function provides a flexible way to construct the quadratic
matrix `W_d` from a desired-gain vector and the observed trait
relationships in a GEBV dataset.

## Usage

``` r
construct_Wd_from_dg(
  gebv_data,
  trait_cols,
  dg,
  lower_is_better = NULL,
  center_traits = TRUE,
  scale_traits = TRUE,
  impute_missing = TRUE,
  method = c("hybrid", "diag_only", "dg_outer", "corr_weighted"),
  base_diag = NULL,
  lambda_diag = 1,
  lambda_outer = 0.5,
  lambda_corr = 1,
  corr_power = 1,
  offdiag_shrink = 1,
  positive_diagonal_only = TRUE,
  normalize = c("trace", "max_abs", "none"),
  debug = TRUE
)
```

## Arguments

- gebv_data:

  A `data.frame` or `data.table` containing trait GEBVs. Each row should
  represent one genotype, and the columns listed in `trait_cols` should
  contain numeric or coercible-to-numeric trait GEBVs.

- trait_cols:

  Character vector giving the names of the trait columns to use when
  constructing `W_d`. These columns must exist in `gebv_data`. Their
  order determines the order used in the returned matrix.

- dg:

  Named numeric desired-gain vector, with names matching `trait_cols`.
  These desired gains are used to build the outer-product and
  correlation-weighted components of `W_d`, and to define default
  diagonal weights when `base_diag = NULL`.

- lower_is_better:

  Optional character vector of trait names for which smaller values are
  favorable. These traits are internally multiplied by `-1` so that all
  traits are oriented in a common favorable direction before
  correlations and matrix components are computed.

- center_traits:

  Logical. If `TRUE`, center the selected GEBV columns by subtracting
  their means before deriving trait relationships. This is useful when
  trait GEBVs are not already centered around zero.

- scale_traits:

  Logical. If `TRUE`, divide each selected trait column by its standard
  deviation after optional centering. This is useful when traits are on
  different scales and you want the correlation structure to be based on
  standardized variables.

- impute_missing:

  Logical. If `TRUE`, missing values are imputed by the mean of the
  corresponding trait column. If `FALSE`, the function stops with an
  error if missing values are present in the selected trait columns.

- method:

  Character string specifying how `W_d` should be constructed. One of:

  - `"hybrid"`

  - `"diag_only"`

  - `"dg_outer"`

  - `"corr_weighted"`

- base_diag:

  Optional named numeric vector giving baseline diagonal weights. Names
  must match `trait_cols`. If `NULL`, the function uses `abs(dg)` as the
  default diagonal weights.

- lambda_diag:

  Numeric multiplier applied to the diagonal component of `W_d`. Larger
  values increase the contribution of squared trait effects.

- lambda_outer:

  Numeric multiplier applied to the desired-gain outer product
  component. Larger values increase the contribution of pairwise trait
  interactions implied by the desired-gain vector alone.

- lambda_corr:

  Numeric multiplier applied to the correlation-weighted component.
  Larger values increase the influence of the empirical trait
  correlation structure on `W_d`.

- corr_power:

  Numeric exponent applied to the absolute trait correlations before
  building the correlation-weighted component. Values greater than `1`
  emphasize strong correlations more heavily; values between `0` and `1`
  reduce contrast among correlation magnitudes.

- offdiag_shrink:

  Numeric shrinkage factor applied to off-diagonal elements of `W_d`.
  Must be between `0` and `1`. A value of `1` leaves off-diagonal
  elements unchanged; smaller values weaken cross-trait interaction
  terms.

- positive_diagonal_only:

  Logical. If `TRUE`, negative diagonal entries are replaced by zero
  after construction. This ensures non-negative diagonal contributions
  in the final matrix.

- normalize:

  Character string controlling normalization of the final matrix. One
  of:

  - `"trace"`: divide `W_d` by the sum of absolute diagonal values;

  - `"max_abs"`: divide `W_d` by its maximum absolute entry;

  - `"none"`: do not normalize.

- debug:

  Logical. If `TRUE`, print progress and diagnostic messages during
  matrix construction.

## Value

A list containing:

- W_d:

  The final symmetric quadratic desired-gain matrix.

- trait_matrix_used:

  A `data.table` containing the transformed GEBV matrix actually used to
  derive trait relationships after optional imputation, direction
  flipping, centering, and scaling.

- cor_matrix:

  The empirical trait correlation matrix computed from the transformed
  GEBV matrix.

- diag_component:

  The diagonal component used in constructing `W_d`.

- outer_component:

  The desired-gain outer-product component.

- corr_component:

  The correlation-weighted component.

- method:

  The construction method used.

## Details

The matrix `W_d` governs the quadratic part of the desired-gain QGSI. It
can be interpreted as controlling how strongly:

- squared trait effects contribute to the index,

- pairwise trait interactions contribute to the index,

- correlation structure among traits influences the quadratic term.

The function supports four strategies for constructing `W_d`:

- `"diag_only"`:

  Builds a purely diagonal matrix, so the quadratic term contains only
  squared trait effects and no cross-trait interactions.

- `"dg_outer"`:

  Uses the outer product of the desired-gain vector, \\d d'\\, to
  construct `W_d`. This introduces pairwise interaction terms based only
  on desired gains.

- `"corr_weighted"`:

  Constructs `W_d` by weighting the desired-gain outer product by the
  empirical trait correlation matrix derived from the transformed GEBV
  data. This allows the quadratic matrix to reflect observed
  relationships among traits.

- `"hybrid"`:

  Combines a diagonal component, an outer-product component, and a
  correlation-weighted component. This is the most flexible option and
  is generally a good default when both squared effects and trait
  interactions are of interest.

The construction workflow is:

1.  convert selected GEBV columns to numeric;

2.  optionally impute missing values;

3.  optionally flip traits in `lower_is_better`;

4.  optionally center and/or scale the GEBV matrix;

5.  compute a trait correlation matrix from the transformed GEBV matrix;

6.  construct the diagonal, outer-product, and correlation-weighted
    components;

7.  combine these components according to `method` and the lambda
    weights;

8.  optionally shrink off-diagonal terms;

9.  enforce symmetry and optional diagonal positivity;

10. optionally normalize the final matrix.

If trait GEBVs are already centered and scaled, you may wish to use:

- `center_traits = FALSE`

- `scale_traits = FALSE`

The returned `W_d` can then be passed directly to
[`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DGQGSI/reference/run_qgsi_desired_gain.md).

## References

Cerón-Rojas JJ, Montesinos-López OA, Montesinos-López A, et al. (2026).
Nonlinear genomic selection index accelerates multi-trait crop
improvement. *Nature Communications*, 17:1991.
doi:10.1038/s41467-026-69890-3

## Examples

``` r
trait_cols <- c("YLD", "MY", "MI", "BL", "NBL", "VHB")

dg <- c(
  YLD = 1.5,
  MY  = 0.5,
  MI  = 0.5,
  BL  = 1.0,
  NBL = 1.0,
  VHB = 1.0
)

ext <- system.file("extdata", package = "DGQGSI")
gebv <- data.table::fread(file.path(ext, "example_gebv.csv"))

Wobj <- construct_Wd_from_dg(
  gebv_data = gebv[, c("GenoID", trait_cols), with = FALSE],
  trait_cols = trait_cols,
  dg = dg,
  lower_is_better = c("BL", "NBL", "VHB"),
  center_traits = FALSE,
  scale_traits = FALSE,
  method = "hybrid",
  lambda_diag = 1.0,
  lambda_outer = 0.5,
  lambda_corr = 1.0,
  corr_power = 1,
  offdiag_shrink = 0.5,
  normalize = "trace",
  debug = FALSE
)

Wobj$W_d
#>            YLD          MY          MI         BL        NBL         VHB
#> YLD 0.34513274 0.010383022 0.013456554 0.03043078 0.03265108 0.027745332
#> MY  0.01038302 0.061946903 0.006802449 0.01108874 0.01560836 0.005314580
#> MI  0.01345655 0.006802449 0.061946903 0.01175799 0.01471109 0.006929337
#> BL  0.03043078 0.011088737 0.011757993 0.17699115 0.01962414 0.011813868
#> NBL 0.03265108 0.015608363 0.014711093 0.01962414 0.17699115 0.011482631
#> VHB 0.02774533 0.005314580 0.006929337 0.01181387 0.01148263 0.176991150
Wobj$cor_matrix
#>              YLD         MY           MI          BL         NBL         VHB
#> YLD  1.000000000 -0.1089062  0.006863548  0.07311311  0.11492869  0.02253708
#> MY  -0.108906175  1.0000000  0.268676726  0.12651361  0.38187249 -0.19972622
#> MI   0.006863548  0.2686767  1.000000000  0.16432659  0.33117673 -0.10849245
#> BL   0.073113110  0.1265136  0.164326593  1.00000000  0.05438195 -0.16625822
#> NBL  0.114928691  0.3818725  0.331176733  0.05438195  1.00000000 -0.17561568
#> VHB  0.022537083 -0.1997262 -0.108492447 -0.16625822 -0.17561568  1.00000000
```
