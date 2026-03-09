# Merge desired-gain index and desired-gain QGSI outputs into one comparison table

Merges the outputs of
[`run_desired_gain_index_Joukhadar2024()`](https://FAkohoue.github.io/DGQGSI/reference/run_desired_gain_index_Joukhadar2024.md)
and
[`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DGQGSI/reference/run_qgsi_desired_gain.md)
into one genotype-level comparison table by ID.

This helper function is intended for situations where both the
desired-gain selection index workflow and the desired-gain QGSI workflow
have been run on the same or overlapping sets of genotypes, and the user
wants to compare:

- genotype scores from both methods,

- genotype ranks from both methods,

- selected versus non-selected status from the DGSI workflow,

- agreement or disagreement between ranking systems.

## Usage

``` r
compare_dg_and_qgsi(
  dg_result,
  qgsi_result,
  id_col = "GenoID",
  dg_cols = c(id_col, "SelectionIndex", "Selected"),
  qgsi_cols = c(id_col, "LinearDGPart", "QuadraticDGPart", "QGSI_DG", "Rank_QGSI_DG"),
  include_metadata = TRUE,
  compute_rank_differences = TRUE,
  add_correlation_summary = TRUE,
  sort_by = c("DG_rank", "QGSI_rank", "DG", "QGSI", "none"),
  debug = TRUE
)
```

## Arguments

- dg_result:

  A result list returned by
  [`run_desired_gain_index_Joukhadar2024()`](https://FAkohoue.github.io/DGQGSI/reference/run_desired_gain_index_Joukhadar2024.md).
  The function expects at least a `ranked_geno` component containing the
  DGSI output table.

- qgsi_result:

  A result list returned by
  [`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DGQGSI/reference/run_qgsi_desired_gain.md).
  The function expects at least a `ranked_geno` component containing the
  QGSI output table.

- id_col:

  Character string naming the genotype identifier column used to merge
  DGSI and QGSI outputs. Default is `"GenoID"`. This column must be
  present in both `dg_result$ranked_geno` and `qgsi_result$ranked_geno`.

- dg_cols:

  Character vector giving the columns to extract from
  `dg_result$ranked_geno`. By default, this includes the genotype
  identifier, the DG selection index, and the selected/non-selected
  flag. Additional metadata columns such as family or population can
  also be requested here.

- qgsi_cols:

  Character vector giving the columns to extract from
  `qgsi_result$ranked_geno`. By default, this includes the genotype
  identifier, the QGSI linear and quadratic components, the final QGSI
  score, and the QGSI rank.

- include_metadata:

  Logical retained for API consistency and future use. At present,
  metadata are included whenever they are explicitly requested in
  `dg_cols` or `qgsi_cols`. Default is `TRUE`.

- compute_rank_differences:

  Logical. If `TRUE`, the function computes diagnostic columns comparing
  DGSI and QGSI ranks, including:

  - `RankDiff_DG_minus_QGSI`

  - `AbsRankDiff_DG_vs_QGSI`

  - `SignSame_DG_QGSI` when both score columns are available

- add_correlation_summary:

  Logical. If `TRUE`, the function computes a summary table containing
  Pearson and Spearman correlations between:

  - DG selection index and QGSI score

  - DG rank and QGSI rank

  whenever the required columns are available.

- sort_by:

  Character string specifying how the merged output table should be
  ordered. One of:

  - `"DG_rank"`: sort by DG rank ascending

  - `"QGSI_rank"`: sort by QGSI rank ascending

  - `"DG"`: sort by DG score descending

  - `"QGSI"`: sort by QGSI score descending

  - `"none"`: do not apply additional sorting after merging

- debug:

  Logical. If `TRUE`, the function prints progress and diagnostic
  messages, including missing requested columns, computed correlation
  summaries, and the final sorting rule applied.

## Value

A list containing:

- comparison_table:

  A `data.table` containing the merged genotype-level outputs from DGSI
  and QGSI, optionally augmented with rank-difference diagnostics.

- correlation_summary:

  A `data.table` containing Pearson and Spearman correlations between
  DGSI and QGSI scores and/or ranks when
  `add_correlation_summary = TRUE`; otherwise `NULL`.

## Details

The function:

1.  extracts selected columns from `dg_result$ranked_geno`;

2.  extracts selected columns from `qgsi_result$ranked_geno`;

3.  renames overlapping DG columns for clarity;

4.  merges both tables by `id_col`;

5.  optionally computes rank-difference diagnostics;

6.  optionally computes correlation summaries between DGSI and QGSI
    scores and/or ranks;

7.  optionally sorts the final merged table according to the chosen
    ranking criterion.

This is especially useful for downstream interpretation, because it
allows side-by-side comparison of:

- `DG_SelectionIndex`

- `QGSI_DG`

- their respective ranks

- and the magnitude of rank disagreement

If some requested columns are not found in either result object, they
are skipped with a debug message rather than causing immediate failure,
as long as the identifier column is available.

## See also

[`run_desired_gain_index_Joukhadar2024()`](https://FAkohoue.github.io/DGQGSI/reference/run_desired_gain_index_Joukhadar2024.md),
[`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DGQGSI/reference/run_qgsi_desired_gain.md),
[`run_dgsi_qgsi_pipeline()`](https://FAkohoue.github.io/DGQGSI/reference/run_dgsi_qgsi_pipeline.md)

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
gebv  <- data.table::fread(file.path(ext, "example_gebv.csv"))

dg_res <- run_desired_gain_index_Joukhadar2024(
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

qgsi_res <- run_qgsi_desired_gain(
  init_data = pheno[, .(GenoID, Family)],
  gebv_data = gebv[, c("GenoID", trait_cols), with = FALSE],
  trait_cols = trait_cols,
  dg = dg,
  lower_is_better = c("BL", "NBL", "VHB"),
  center_traits = FALSE,
  scale_traits = FALSE,
  debug = FALSE
)

merged <- compare_dg_and_qgsi(
  dg_result = dg_res,
  qgsi_result = qgsi_res,
  id_col = "GenoID",
  compute_rank_differences = TRUE,
  add_correlation_summary = TRUE,
  sort_by = "DG_rank",
  debug = FALSE
)

head(merged$comparison_table)
#>    GenoID DG_SelectionIndex DG_Selected LinearDGPart QuadraticDGPart  QGSI_DG
#>    <char>             <num>      <lgcl>        <num>           <num>    <num>
#> 1:   G013          9.095305       FALSE     3.043771       2.4040859 5.447857
#> 2:   G016          7.375457       FALSE     2.934136       2.9481011 5.882237
#> 3:   G029          7.332733       FALSE     3.245327       3.4324322 6.677759
#> 4:   G001          6.995980        TRUE     2.581310       1.2789864 3.860296
#> 5:   G007          6.072997       FALSE     1.752643       0.8931555 2.645798
#> 6:   G025          5.851942       FALSE     1.849823       4.6185705 6.468394
#>    Rank_QGSI_DG SignSame_DG_QGSI
#>           <num>           <lgcl>
#> 1:            4             TRUE
#> 2:            3             TRUE
#> 3:            1             TRUE
#> 4:            5             TRUE
#> 5:           11             TRUE
#> 6:            2             TRUE
merged$correlation_summary
#>                      Comparison   Pearson  Spearman     N
#>                          <char>     <num>     <num> <int>
#> 1: DG_SelectionIndex vs QGSI_DG 0.8510171 0.8874305    30
```
