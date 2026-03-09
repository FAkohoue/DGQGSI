# Run desired-gain QGSI using trait GEBVs

Computes a desired-gain adaptation of the quadratic genomic selection
index.

## Usage

``` r
run_qgsi_desired_gain(init_data, gebv_data, trait_cols, id_col = "GenoID",
  dg, W_d = NULL, quadratic_diag_weights = NULL, lower_is_better = NULL,
  center_traits = TRUE, scale_traits = TRUE, impute_missing = TRUE,
  return_components = TRUE, debug = TRUE)
```

## Arguments

- init_data:

  Metadata table.

- gebv_data:

  Trait GEBV table.

- trait_cols:

  Trait columns.

- id_col:

  Identifier column.

- dg:

  Desired-gain vector.

- W_d:

  Quadratic matrix.

- quadratic_diag_weights:

  Diagonal weights if W_d is NULL.

- lower_is_better:

  Traits to flip.

- center_traits:

  Center GEBVs.

- scale_traits:

  Scale GEBVs.

- impute_missing:

  Impute missing GEBVs.

- return_components:

  Return linear and quadratic parts.

- debug:

  Print debug messages.

## Value

A list containing the transformed trait space, ranked genotypes, and
component summaries.
