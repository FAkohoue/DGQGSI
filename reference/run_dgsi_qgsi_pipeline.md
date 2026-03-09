# Run a full DG and QGSI pipeline

High-level package wrapper for one-stop execution of the DG and
desired-gain QGSI workflows.

## Usage

``` r
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
```

## Arguments

- mode:

  Run DG only, QGSI only, or both.

- init_data:

  Metadata table.

- trait_cols:

  Trait columns.

- dg:

  Desired-gain vector.

- id_col:

  Identifier column.

- cand_data:

  Candidate trait table for DG.

- ref_data:

  Reference table for DG.

- G:

  Genetic covariance matrix for DG.

- select_mode:

  Selection mode for DG.

- n_select:

  Number selected in top-n mode.

- trait_min_sd:

  Thresholds for DG.

- fallback_to_top_n:

  Fallback behavior.

- n_iter:

  DG iterations.

- n_rep:

  DG replicates.

- sd_scale:

  DG perturbation scale.

- seed:

  Random seed.

- ridge_P:

  DG ridge for P.

- ridge_M:

  DG ridge for M.

- dg_scale_traits:

  Scale DG inputs.

- gebv_data:

  Trait GEBV table for QGSI.

- W_d:

  Optional quadratic matrix.

- quadratic_diag_weights:

  Optional diagonal weights.

- qgsi_center_traits:

  Center GEBVs.

- qgsi_scale_traits:

  Scale GEBVs.

- qgsi_impute_missing:

  Impute missing GEBVs.

- W_method:

  Automatic W_d method.

- W_base_diag:

  Optional base diagonal.

- W_lambda_diag:

  Diagonal weight.

- W_lambda_outer:

  Outer-product weight.

- W_lambda_corr:

  Correlation-weighted component weight.

- W_corr_power:

  Correlation exponent.

- W_offdiag_shrink:

  Off-diagonal shrinkage.

- W_positive_diagonal_only:

  Force non-negative diagonal.

- W_normalize:

  Normalization method.

- lower_is_better:

  Traits to flip.

- merge_outputs:

  Merge outputs when mode is both.

- compare_sort_by:

  Sort order for merged output.

- debug:

  Print debug messages.

## Value

A list containing dg_result, W_d_result, qgsi_result, and optionally
comparison_result.
