# Compare DG and QGSI outputs

Merges DG and QGSI outputs into a single comparison table by genotype
ID.

## Usage

``` r
compare_dg_and_qgsi(dg_result, qgsi_result, id_col = "GenoID",
  dg_cols = c(id_col, "SelectionIndex", "Selected"),
  qgsi_cols = c(id_col, "LinearDGPart", "QuadraticDGPart", "QGSI_DG", "Rank_QGSI_DG"),
  include_metadata = TRUE, compute_rank_differences = TRUE,
  add_correlation_summary = TRUE, sort_by = c("DG_rank", "QGSI_rank", "DG", "QGSI", "none"),
  debug = TRUE)
```

## Arguments

- dg_result:

  Output of run_desired_gain_index_Joukhadar2024().

- qgsi_result:

  Output of run_qgsi_desired_gain().

- id_col:

  Identifier column.

- dg_cols:

  Columns to extract from DG result.

- qgsi_cols:

  Columns to extract from QGSI result.

- include_metadata:

  Retained for API consistency.

- compute_rank_differences:

  Compute rank differences.

- add_correlation_summary:

  Compute score/rank correlations.

- sort_by:

  Final sort order.

- debug:

  Print debug messages.

## Value

A list with comparison_table and optional correlation_summary.
