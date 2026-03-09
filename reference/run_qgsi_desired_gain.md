# Desired-gain QGSI using trait GEBVs

Builds a genomic quadratic selection index for all genotypes using the
QGSI structure proposed for nonlinear genomic selection, but replacing
the linear economic-weight vector with a breeder-defined desired-gain
vector.

The implemented index is: \$\$QGSI\\DG_i = d' \hat{\gamma}\_i +
\hat{\gamma}\_i' W_d \hat{\gamma}\_i\$\$

where:

- \\d\\ is the desired-gain vector,

- \\\hat{\gamma}\_i\\ is the vector of trait GEBVs for genotype \\i\\,

- \\W_d\\ is a quadratic desired-gain matrix controlling trait-by-trait
  interaction contributions.

## Usage

``` r
run_qgsi_desired_gain(
  init_data,
  gebv_data,
  trait_cols,
  id_col = "GenoID",
  dg,
  W_d = NULL,
  quadratic_diag_weights = NULL,
  lower_is_better = NULL,
  center_traits = TRUE,
  scale_traits = TRUE,
  impute_missing = TRUE,
  return_components = TRUE,
  debug = TRUE
)
```

## Arguments

- init_data:

  A `data.frame` or `data.table` containing genotype identifiers and any
  metadata to preserve in the final output. Typical columns may include
  `GenoID`, family, population, or other descriptors. The returned table
  appends QGSI results to this object.

- gebv_data:

  A `data.frame` or `data.table` containing genotype IDs and trait
  GEBVs. It must contain `id_col` and all columns listed in
  `trait_cols`. Each row should represent one genotype and each trait
  column should be numeric or coercible to numeric.

- trait_cols:

  Character vector giving the trait GEBV column names to use in the
  index. These names must be present in `gebv_data`. The order of
  `trait_cols` determines the order used for `dg`, `W_d`, and any
  diagonal weights.

- id_col:

  Character string naming the genotype identifier column used to align
  `init_data` and `gebv_data`. Default is `"GenoID"`. This column must
  be present in both input tables.

- dg:

  Named numeric desired-gain vector, with names matching `trait_cols`.
  This replaces the classical economic-weight vector in the linear part
  of the index. For best interpretability, `dg` should be on the same
  scale as the transformed GEBV matrix used by the function.

- W_d:

  Optional full symmetric quadratic desired-gain matrix of dimension
  `length(trait_cols) x length(trait_cols)`. This matrix controls the
  quadratic and interaction contributions among traits. If `NULL`, a
  diagonal matrix is created from `quadratic_diag_weights`.

- quadratic_diag_weights:

  Optional named numeric vector of diagonal weights used only when
  `W_d = NULL`. Names must match `trait_cols`. Larger values place
  stronger emphasis on the quadratic contribution of the corresponding
  trait.

- lower_is_better:

  Optional character vector naming traits for which smaller values are
  favorable. These traits are internally multiplied by `-1` so that all
  traits are aligned in a common favorable direction before building the
  index.

- center_traits:

  Logical. If `TRUE`, the GEBV matrix is centered trait by trait using
  the column means before the index is computed. This is useful when
  GEBVs are not already centered around zero.

- scale_traits:

  Logical. If `TRUE`, each trait column is divided by its standard
  deviation after optional centering. This is useful when GEBVs are on
  different scales or when you want a standardized trait space.

- impute_missing:

  Logical. If `TRUE`, missing trait values are imputed by the mean of
  the corresponding trait column. If `FALSE`, the function stops with an
  error when missing values are present.

- return_components:

  Logical. If `TRUE`, the returned ranked table includes the separate
  linear and quadratic components in addition to the final QGSI score.
  If `FALSE`, only the final QGSI score and its rank are retained.

- debug:

  Logical. If `TRUE`, the function prints progress and diagnostic
  messages during execution. This is useful for tracing data preparation
  and matrix construction steps.

## Value

A list containing:

- dg:

  The desired-gain vector used in the linear component.

- W_d:

  The quadratic desired-gain matrix used in the quadratic component.

- trait_space_used:

  A `data.table` containing the transformed GEBV matrix actually used in
  the index calculation after optional flipping, centering, and scaling.

- ranked_geno:

  A `data.table` containing `init_data` augmented with `QGSI_DG`,
  `Rank_QGSI_DG`, and optionally `LinearDGPart`, `QuadraticDGPart`,
  `Rank_LinearDGPart`, and `Rank_QuadraticDGPart`, sorted by decreasing
  `QGSI_DG`.

- component_summary:

  A `data.table` summarizing the mean, standard deviation, minimum, and
  maximum of the computed index components.

## Details

This function is intended for genomic selection workflows where trait
GEBVs are already available for each genotype.

The function:

1.  aligns `init_data` and `gebv_data` by `id_col`;

2.  converts all requested trait columns to numeric;

3.  optionally imputes missing values by trait means;

4.  optionally flips traits listed in `lower_is_better` so that all
    traits are expressed in a common favorable direction;

5.  optionally centers and/or scales the trait GEBVs;

6.  computes a linear desired-gain component;

7.  computes a quadratic interaction component using `W_d`;

8.  combines both parts into a final desired-gain QGSI score;

9.  ranks all genotypes by the resulting index.

If traits are already standardized and `dg` is already expressed in SD
units, the recommended settings are:

- `center_traits = FALSE`

- `scale_traits = FALSE`

because the trait GEBVs and desired gains are already on the same scale.

If `W_d` is not supplied, the function builds a diagonal matrix from
`quadratic_diag_weights`. In that case, the quadratic contribution
captures only squared trait effects, without off-diagonal trait
interactions.

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

qdiag <- c(
  YLD = 2.0,
  MY  = 1.0,
  MI  = 1.0,
  BL  = 1.5,
  NBL = 1.5,
  VHB = 1.5
)

ext <- system.file("extdata", package = "DGQGSI")
pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
gebv  <- data.table::fread(file.path(ext, "example_gebv.csv"))

res <- run_qgsi_desired_gain(
  init_data = pheno[, .(GenoID, Family)],
  gebv_data = gebv[, c("GenoID", trait_cols), with = FALSE],
  trait_cols = trait_cols,
  dg = dg,
  quadratic_diag_weights = qdiag,
  lower_is_better = c("BL", "NBL", "VHB"),
  center_traits = FALSE,
  scale_traits = FALSE,
  impute_missing = TRUE,
  return_components = TRUE,
  debug = FALSE
)

head(res$ranked_geno)
#>    GenoID Family LinearDGPart QuadraticDGPart   QGSI_DG Rank_QGSI_DG
#>    <char> <char>        <num>           <num>     <num>        <num>
#> 1:   G025     F5     1.849823        8.276079 10.125902            1
#> 2:   G029     F4     3.245327        5.339413  8.584740            2
#> 3:   G016     F1     2.934136        4.459293  7.393430            3
#> 4:   G013     F3     3.043771        4.127950  7.171721            4
#> 5:   G001     F1     2.581310        2.076364  4.657674            5
#> 6:   G015     F5     1.281385        3.313659  4.595043            6
#>    Rank_LinearDGPart Rank_QuadraticDGPart
#>                <num>                <num>
#> 1:                 6                    1
#> 2:                 1                    2
#> 3:                 3                    4
#> 4:                 2                    6
#> 5:                 4                   17
#> 6:                 9                    9
res$component_summary
#>          Component      Mean       SD        Min       Max
#>             <char>     <num>    <num>      <num>     <num>
#> 1:    LinearDGPart 0.5528237 1.448067 -2.6519349  3.245327
#> 2: QuadraticDGPart 2.6898094 1.628205  0.4977952  8.276079
#> 3:         QGSI_DG 3.2426331 2.361725  0.2290951 10.125902
```
