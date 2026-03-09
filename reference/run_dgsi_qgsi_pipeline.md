# Run a full desired-gain index and desired-gain QGSI pipeline

High-level wrapper that runs one or both of the desired-gain selection
index optimization workflow and the desired-gain quadratic genomic
selection index (QGSI) workflow, with optional merged comparison output.

This function orchestrates the full DGQGSI pipeline, allowing users to:

- compute a desired-gain selection index following Joukhadar et al.
  (2024),

- compute a desired-gain quadratic genomic selection index (QGSI)
  following Cerón-Rojas et al. (2026),

- automatically construct the quadratic desired-gain matrix `W_d` when
  needed,

- and generate a merged comparison table between the two methods.

The function acts as a convenient wrapper that coordinates multiple
package functions, ensuring consistent input alignment and reducing
repetitive code when running complete selection-index analyses.

## Usage

``` r
run_dgsi_qgsi_pipeline(
  mode = c("both", "dg", "qgsi"),
  init_data,
  trait_cols,
  dg,
  id_col = "GenoID",
  cand_data = NULL,
  ref_data = NULL,
  G = NULL,
  select_mode = c("top_n", "trait_thresholds"),
  n_select = 100,
  trait_min_sd = NULL,
  fallback_to_top_n = TRUE,
  n_iter = 1000,
  n_rep = 20,
  sd_scale = 1,
  seed = 42,
  ridge_P = 1e-06,
  ridge_M = 1e-06,
  dg_scale_traits = FALSE,
  gebv_data = NULL,
  W_d = NULL,
  quadratic_diag_weights = NULL,
  qgsi_center_traits = FALSE,
  qgsi_scale_traits = FALSE,
  qgsi_impute_missing = TRUE,
  W_method = c("hybrid", "diag_only", "dg_outer", "corr_weighted"),
  W_base_diag = NULL,
  W_lambda_diag = 1,
  W_lambda_outer = 0.5,
  W_lambda_corr = 1,
  W_corr_power = 1,
  W_offdiag_shrink = 1,
  W_positive_diagonal_only = TRUE,
  W_normalize = c("trace", "max_abs", "none"),
  lower_is_better = NULL,
  merge_outputs = TRUE,
  compare_sort_by = c("DG_rank", "QGSI_rank", "DG", "QGSI", "none"),
  debug = TRUE
)
```

## Arguments

- mode:

  Character string specifying which workflow to run. One of `"both"`,
  `"dg"`, or `"qgsi"`.

- init_data:

  A `data.frame` or `data.table` containing genotype identifiers and
  metadata (for example `GenoID`, family, population). This table is
  preserved and appended with index results.

- trait_cols:

  Character vector of trait names used in both DG and QGSI calculations.
  These columns must exist in the candidate trait data (`cand_data`)
  and/or the GEBV data (`gebv_data`).

- dg:

  Named numeric desired-gain vector whose names correspond to
  `trait_cols`. This vector defines the breeding objectives and is used
  in both DGSI and QGSI calculations.

- id_col:

  Character genotype identifier column used to align all input tables.
  Default is `"GenoID"`.

- cand_data:

  Optional `data.frame` or `data.table` containing candidate trait
  values for the desired-gain index workflow.

- ref_data:

  Optional reference dataset used for scaling and covariance estimation
  in the desired-gain index workflow.

- G:

  Optional genetic covariance matrix used in the desired-gain index
  optimization. If `NULL`, it is estimated from the reference data.

- select_mode:

  Selection strategy used in the DG workflow. One of `"top_n"` or
  `"trait_thresholds"`.

- n_select:

  Integer specifying the number of selected genotypes when
  `select_mode = "top_n"`.

- trait_min_sd:

  Named numeric vector specifying minimum selection thresholds when
  `select_mode = "trait_thresholds"`. Values are expressed in standard
  deviation units.

- fallback_to_top_n:

  Logical. If `TRUE`, the function falls back to top-N selection when
  threshold selection returns zero genotypes.

- n_iter:

  Integer number of optimization iterations used in the desired-gain
  index search algorithm.

- n_rep:

  Integer number of optimization replicates.

- sd_scale:

  Numeric perturbation scale applied to candidate desired-gain vectors
  during optimization.

- seed:

  Integer random seed controlling reproducibility of the optimization
  procedure.

- ridge_P:

  Numeric ridge penalty added to the phenotypic correlation matrix used
  in the desired-gain index calculation.

- ridge_M:

  Numeric ridge penalty added to the inner matrix used during
  coefficient computation.

- dg_scale_traits:

  Logical. If `TRUE`, traits are centered and scaled before computing
  the desired-gain index.

- gebv_data:

  Optional `data.frame` or `data.table` containing genotype-level GEBVs
  used for QGSI computation.

- W_d:

  Optional quadratic desired-gain matrix. If `NULL`, the matrix is
  automatically generated using
  [`construct_Wd_from_dg()`](https://FAkohoue.github.io/DGQGSI/reference/construct_Wd_from_dg.md).

- quadratic_diag_weights:

  Optional named numeric vector used to construct a diagonal `W_d`
  matrix when no full matrix is provided.

- qgsi_center_traits:

  Logical. If `TRUE`, the GEBV matrix is centered before computing QGSI.

- qgsi_scale_traits:

  Logical. If `TRUE`, the GEBV matrix is scaled before computing QGSI.

- qgsi_impute_missing:

  Logical. If `TRUE`, missing GEBV values are imputed by trait means
  before computing QGSI.

- W_method:

  Method used to automatically construct `W_d`. Passed to
  [`construct_Wd_from_dg()`](https://FAkohoue.github.io/DGQGSI/reference/construct_Wd_from_dg.md).
  Options include `"hybrid"`, `"diag_only"`, `"dg_outer"`, and
  `"corr_weighted"`.

- W_base_diag:

  Optional baseline diagonal weights for automatic `W_d` construction.

- W_lambda_diag:

  Weight applied to the diagonal component when building `W_d`.

- W_lambda_outer:

  Weight applied to the desired-gain outer-product component of `W_d`.

- W_lambda_corr:

  Weight applied to the correlation-weighted component of `W_d`.

- W_corr_power:

  Exponent applied to trait correlations when constructing the
  correlation-weighted component.

- W_offdiag_shrink:

  Numeric shrinkage factor applied to off-diagonal elements of `W_d`.

- W_positive_diagonal_only:

  Logical. If `TRUE`, negative diagonal entries in `W_d` are set to
  zero.

- W_normalize:

  Normalization strategy applied to `W_d`. One of `"trace"`,
  `"max_abs"`, or `"none"`.

- lower_is_better:

  Optional character vector listing traits for which smaller values are
  favorable. These traits are internally flipped so that all objectives
  align in the same direction.

- merge_outputs:

  Logical. If `TRUE` and `mode = "both"`, the DG and QGSI outputs are
  merged using
  [`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DGQGSI/reference/compare_dg_and_qgsi.md).

- compare_sort_by:

  Sorting rule passed to
  [`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DGQGSI/reference/compare_dg_and_qgsi.md).
  Determines how the merged comparison table is ordered.

- debug:

  Logical. If `TRUE`, the function prints progress and diagnostic
  messages during execution.

## Value

A list that may contain:

- dg_result:

  Output of
  [`run_desired_gain_index_Joukhadar2024()`](https://FAkohoue.github.io/DGQGSI/reference/run_desired_gain_index_Joukhadar2024.md).

- W_d_result:

  Output of
  [`construct_Wd_from_dg()`](https://FAkohoue.github.io/DGQGSI/reference/construct_Wd_from_dg.md)
  when the quadratic matrix is constructed automatically.

- qgsi_result:

  Output of
  [`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DGQGSI/reference/run_qgsi_desired_gain.md).

- comparison_result:

  Output of
  [`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DGQGSI/reference/compare_dg_and_qgsi.md)
  when both workflows are run and `merge_outputs = TRUE`.

## Details

The pipeline can operate in three modes controlled by `mode`:

- `"dg"`:

  Runs only the desired-gain selection index workflow using
  [`run_desired_gain_index_Joukhadar2024()`](https://FAkohoue.github.io/DGQGSI/reference/run_desired_gain_index_Joukhadar2024.md).

- `"qgsi"`:

  Runs only the desired-gain QGSI workflow using
  [`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DGQGSI/reference/run_qgsi_desired_gain.md).
  If `W_d` is not provided, the matrix is automatically constructed
  using
  [`construct_Wd_from_dg()`](https://FAkohoue.github.io/DGQGSI/reference/construct_Wd_from_dg.md).

- `"both"`:

  Runs both workflows sequentially. If `merge_outputs = TRUE`, the
  results are merged using
  [`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DGQGSI/reference/compare_dg_and_qgsi.md)
  to produce a genotype-level comparison table.

The typical pipeline sequence when `mode = "both"` is:

1.  run the desired-gain index optimization;

2.  optionally construct a quadratic desired-gain matrix `W_d`;

3.  compute the QGSI scores using GEBVs;

4.  optionally merge both outputs into a comparison table.

This wrapper is particularly useful in genomic selection pipelines where
breeders want to compare:

- classical linear desired-gain indices,

- nonlinear quadratic genomic selection indices,

- and their resulting genotype rankings.

## References

Joukhadar R et al. (2024). Optimising desired gain indices to maximise
selection response. Frontiers in Plant Science 15:1337388.

Cerón-Rojas JJ et al. (2026). Nonlinear genomic selection index
accelerates multi-trait crop improvement. Nature Communications 17:1991.

## See also

[`run_desired_gain_index_Joukhadar2024()`](https://FAkohoue.github.io/DGQGSI/reference/run_desired_gain_index_Joukhadar2024.md),
[`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DGQGSI/reference/run_qgsi_desired_gain.md),
[`construct_Wd_from_dg()`](https://FAkohoue.github.io/DGQGSI/reference/construct_Wd_from_dg.md),
[`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DGQGSI/reference/compare_dg_and_qgsi.md)

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
pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
gebv  <- data.table::fread(file.path(ext, "example_gebv.csv"))

res <- run_dgsi_qgsi_pipeline(
  mode = "both",
  init_data = pheno[, .(GenoID, Family)],
  cand_data = pheno[, c("GenoID", trait_cols), with = FALSE],
  gebv_data = gebv[, c("GenoID", trait_cols), with = FALSE],
  trait_cols = trait_cols,
  dg = dg,
  lower_is_better = c("BL", "NBL", "VHB"),
  merge_outputs = TRUE,
  debug = FALSE
)

head(res$comparison_result$comparison_table)
#>    GenoID DG_SelectionIndex DG_Selected LinearDGPart QuadraticDGPart  QGSI_DG
#>    <char>             <num>      <lgcl>        <num>           <num>    <num>
#> 1:   G013          9.095305        TRUE     3.043771       0.8199342 3.863705
#> 2:   G016          7.375457        TRUE     2.934136       0.7698778 3.704014
#> 3:   G029          7.332733        TRUE     3.245327       0.9101123 4.155439
#> 4:   G001          6.995980        TRUE     2.581310       0.4790381 3.060348
#> 5:   G007          6.072997        TRUE     1.752643       0.2413784 1.994021
#> 6:   G025          5.851942        TRUE     1.849823       1.1382887 2.988112
#>    Rank_QGSI_DG SignSame_DG_QGSI
#>           <num>           <lgcl>
#> 1:            2             TRUE
#> 2:            3             TRUE
#> 3:            1             TRUE
#> 4:            4             TRUE
#> 5:            7             TRUE
#> 6:            5             TRUE
```
