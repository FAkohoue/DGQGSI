# Run a desired-gain selection index

Optimizes desired-gain index coefficients and returns a ranked genotype
table plus selected subsets.

## Usage

``` r
run_desired_gain_index_Joukhadar2024(init_data, cand_data, trait_cols, ref_data = NULL,
  id_col = "GenoID", scale_traits = FALSE, lower_is_better = NULL,
  G = NULL, dg, select_mode = c("top_n", "trait_thresholds"),
  n_select = 100, trait_min_sd = NULL, fallback_to_top_n = TRUE,
  n_iter = 1000, n_rep = 20, sd_scale = 1, seed = 42,
  ridge_P = 1e-06, ridge_M = 1e-06, debug = TRUE, return_all_reps = TRUE)
```

## Arguments

- init_data:

  Metadata table.

- cand_data:

  Candidate trait table.

- trait_cols:

  Trait columns.

- ref_data:

  Reference table.

- id_col:

  Identifier column.

- scale_traits:

  Scale traits or not.

- lower_is_better:

  Traits to flip.

- G:

  Genetic covariance matrix.

- dg:

  Desired-gain vector.

- select_mode:

  Selection mode.

- n_select:

  Number selected in top-n mode.

- trait_min_sd:

  Thresholds in threshold mode.

- fallback_to_top_n:

  Fallback behavior.

- n_iter:

  Iterations.

- n_rep:

  Replicates.

- sd_scale:

  Perturbation scale.

- seed:

  Random seed.

- ridge_P:

  Ridge penalty for P.

- ridge_M:

  Ridge penalty for M.

- debug:

  Print debug messages.

- return_all_reps:

  Return all replicates.

## Value

A list with optimized coefficients, realized gains, ranked output, and
optional replicate details.
