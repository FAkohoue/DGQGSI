# Run a desired-gain selection index with optional trait scaling and threshold-based or top-N selection

Implements a desired-gain selection index optimization strategy inspired
by Joukhadar et al. (2024). The function searches for index coefficients
that produce realized multi-trait gains close to a user-defined
desired-gain vector and then ranks all genotypes using the best
replicate.

## Usage

``` r
run_desired_gain_index_Joukhadar2024(
  init_data,
  cand_data,
  trait_cols,
  ref_data = NULL,
  id_col = "GenoID",
  scale_traits = FALSE,
  lower_is_better = NULL,
  G = NULL,
  dg,
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
  debug = TRUE,
  return_all_reps = TRUE
)
```

## Arguments

- init_data:

  A `data.frame` or `data.table` containing genotype identifiers and any
  metadata you want to preserve in the final output. Typical columns may
  include `GenoID`, family, group, or other descriptive information. The
  function returns this table augmented with the final `SelectionIndex`
  and `Selected` columns.

- cand_data:

  A `data.frame` or `data.table` containing the candidate trait values
  used to build the desired-gain index. It must contain `id_col` and all
  columns listed in `trait_cols`. Each row represents one genotype and
  each trait column should be numeric or coercible to numeric.

- trait_cols:

  Character vector giving the names of the trait columns to include in
  the index. These names must be present in both `cand_data` and
  `ref_data` (or in `cand_data` alone if `ref_data = NULL`). The order
  of `trait_cols` defines the order used for `dg`, `trait_min_sd`, `P`,
  and `G`.

- ref_data:

  Optional `data.frame` or `data.table` used as the reference population
  for scaling traits and for estimating the phenotypic correlation
  matrix `P` and, if necessary, the genetic covariance matrix `G`. If
  `NULL`, `cand_data` is used as the reference.

- id_col:

  Character string naming the genotype identifier column used to align
  `init_data` and `cand_data`. Default is `"GenoID"`. The function
  requires that this column exists in both `init_data` and `cand_data`.

- scale_traits:

  Logical. If `TRUE`, the trait columns in both candidate and reference
  data are centered and scaled using the mean and standard deviation of
  `ref_data`. Use this when traits are on different raw scales. If
  `FALSE`, the function assumes traits are already on a comparable
  scale, such as SD units.

- lower_is_better:

  Optional character vector naming traits for which smaller values are
  favorable. These traits are internally multiplied by `-1` so that all
  traits are oriented in a common “higher is better” direction during
  optimization.

- G:

  Optional square genetic covariance matrix among traits. Its dimension
  must be `length(trait_cols) x length(trait_cols)`. If `NULL`, the
  function estimates `G` as `cov(ref_mat)` from the reference trait
  matrix.

- dg:

  Named numeric vector of desired gains, with names matching
  `trait_cols`. Each value represents the desired response for one trait
  in the trait space used by the function. For best interpretability,
  `dg` should be on the same scale as the transformed trait values.

- select_mode:

  Character string specifying how the selected subset is determined. One
  of:

  - `"top_n"`: selects the `n_select` highest-ranked genotypes by the
    final index;

  - `"trait_thresholds"`: selects genotypes meeting minimum trait
    thresholds in the transformed trait space.

- n_select:

  Integer giving the number of genotypes to select when
  `select_mode = "top_n"`. Ignored when
  `select_mode = "trait_thresholds"` unless fallback is triggered.

- trait_min_sd:

  Named numeric vector of minimum acceptable trait values used when
  `select_mode = "trait_thresholds"`. Names must match `trait_cols`.
  Thresholds are applied after direction flipping and optional scaling,
  so they should be specified on that transformed scale.

- fallback_to_top_n:

  Logical. Relevant only when `select_mode = "trait_thresholds"`. If
  `TRUE` and no genotype passes the thresholds, the function falls back
  to selecting the top `n_select` genotypes by the final index. If
  `FALSE`, the function stops with an error.

- n_iter:

  Integer giving the number of optimization iterations performed within
  each replicate. Larger values increase the search effort but also
  increase computation time.

- n_rep:

  Integer giving the number of independent optimization replicates.
  Replicates are useful because the search is stochastic, and different
  replicates may converge to slightly different solutions.

- sd_scale:

  Numeric factor controlling the perturbation scale used when sampling
  candidate desired-gain vectors around the current best solution.
  Larger values increase exploration; smaller values keep the search
  more local around the current optimum.

- seed:

  Integer random seed used to make the stochastic optimization
  reproducible.

- ridge_P:

  Numeric ridge penalty added to the phenotypic correlation matrix `P`
  before inversion. This can improve numerical stability when `P` is
  close to singular.

- ridge_M:

  Numeric ridge penalty added to the inner matrix `M` during the
  coefficient computation. This can improve stability in cases of highly
  correlated traits or poorly conditioned matrices.

- debug:

  Logical. If `TRUE`, the function prints progress and diagnostic
  messages during execution. This is helpful for understanding what the
  function is doing and for troubleshooting.

- return_all_reps:

  Logical. If `TRUE`, the returned object includes the full list of
  replicate-level optimization results. If `FALSE`, only the summary
  corresponding to the best replicate is returned.

## Value

A list containing:

- dg:

  The desired-gain vector used in the analysis.

- trait_min_sd:

  The threshold vector used when `select_mode = "trait_thresholds"`,
  otherwise `NULL`.

- optimized_d:

  The best sampled desired-gain vector from the winning replicate.

- optimized_b:

  The optimized index coefficient vector from the winning replicate.

- realized_g:

  The realized gains achieved by the selected subset in the transformed
  trait space.

- best_q:

  The objective value of the best replicate. Lower values indicate a
  closer match between realized and desired gains.

- q_trace:

  The iteration-wise objective values from the best replicate.

- mean_replicate_cor:

  Mean pairwise correlation among replicate index vectors, used as a
  simple replicate-stability diagnostic.

- best_replicate:

  The index of the replicate with the best objective value.

- ranked_geno:

  A `data.table` containing `init_data` plus the final `SelectionIndex`
  and `Selected` columns, sorted by decreasing `SelectionIndex`.

- selection_summary:

  A `data.table` summarizing mean trait values and realized gains for
  all candidates versus the selected subset.

- selected_geno:

  Subset of `ranked_geno` with `Selected == TRUE`.

- non_selected_geno:

  Subset of `ranked_geno` with `Selected == FALSE`.

- all_reps:

  Optional list of all replicate results when `return_all_reps = TRUE`.

## Details

The workflow:

1.  aligns `init_data` and `cand_data` by `id_col`;

2.  optionally flips traits in `lower_is_better` so all objectives are
    in the same direction;

3.  optionally scales traits using `ref_data`;

4.  computes a phenotypic correlation matrix `P` and a genetic
    covariance matrix `G`;

5.  repeatedly samples candidate desired-gain vectors, computes index
    coefficients, evaluates realized gains on the selected set, and
    keeps the best solution according to a goodness-of-fit criterion.

This function is **selection-oriented**. It does not only return one
index value per genotype, but also determines a selected subset based on
either:

- the highest index values (`select_mode = "top_n"`), or

- user-defined minimum trait thresholds
  (`select_mode = "trait_thresholds"`).

If your trait values are already expressed in comparable standardized
units (for example SD units), set `scale_traits = FALSE`. If raw trait
values are used, `scale_traits = TRUE` is generally more appropriate.

## References

Joukhadar R, Li Y, Thistlethwaite R, Forrest KL, Tibbits JF, Trethowan
R, Hayden MJ (2024). Optimising desired gain indices to maximise
selection response. *Frontiers in Plant Science*, 15:1337388.
doi:10.3389/fpls.2024.1337388

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

trait_min_sd <- c(
  YLD = 0.2,
  MY  = 0.1,
  MI  = 0.1,
  BL  = 0.1,
  NBL = 0.1,
  VHB = 0.1
)

ext <- system.file("extdata", package = "DGQGSI")
pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))

res <- run_desired_gain_index_Joukhadar2024(
  init_data = pheno[, .(GenoID, Family)],
  cand_data = pheno[, c("GenoID", trait_cols), with = FALSE],
  ref_data = pheno[, c("GenoID", trait_cols), with = FALSE],
  trait_cols = trait_cols,
  dg = dg,
  lower_is_better = c("BL", "NBL", "VHB"),
  select_mode = "trait_thresholds",
  trait_min_sd = trait_min_sd,
  scale_traits = FALSE,
  debug = FALSE
)

head(res$ranked_geno)
#>    GenoID Family SelectionIndex Selected
#>    <char> <char>          <num>   <lgcl>
#> 1:   G013     F3       9.095305    FALSE
#> 2:   G016     F1       7.375457    FALSE
#> 3:   G029     F4       7.332733    FALSE
#> 4:   G001     F1       6.995980     TRUE
#> 5:   G007     F2       6.072997    FALSE
#> 6:   G025     F5       5.851942    FALSE
res$selection_summary
#>     Trait    Mean_All Mean_Selected    SD_All SD_Selected Realized_Gain_SD
#>    <char>       <num>         <num>     <num>       <num>            <num>
#> 1:    YLD  0.49133893     1.5881848 1.0801658          NA        1.0154421
#> 2:     MY  0.15544154     1.1595565 0.8077812          NA        1.2430532
#> 3:     MI  0.21967789     0.1464351 0.6253735          NA       -0.1171185
#> 4:     BL  0.14404814     0.8116349 0.9487245          NA        0.7036677
#> 5:    NBL -0.03033278     0.8645213 0.6348214          NA        1.4096154
#> 6:    VHB -0.02779384     0.1217989 0.9537685          NA        0.1568439
```
