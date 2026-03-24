# DGQGSI Introduction

## Overview

`DGQGSI` provides tools for **multi-trait selection in plant breeding
programs** using desired-gain–based selection indices. It integrates two
complementary workflows:

- **DGSI** — Desired-Gain Selection Index Optimization: estimates index
  weights whose realized genetic gains closely match breeder-defined
  targets (Joukhadar et al., 2024).
- **QGSI** — Desired-Gain Quadratic Genomic Selection Index: a nonlinear
  index that scores genotypes using both linear and quadratic GEBV
  components, capturing trait interactions that a linear index would
  miss (Cerón-Rojas et al., 2026).

The key advantage of both methods is that breeders specify **desired
genetic gains** directly — such as *“increase yield by 1.5 units while
reducing disease score by 1 unit”* — rather than economic weights, which
are rarely known with precision in practice.

------------------------------------------------------------------------

## Methodological background

### DGSI

The DGSI finds a weight vector **b** such that the linear index

$$I = \mathbf{b}^{\top}\mathbf{X}$$

produces realized genetic gains **g** as close as possible to the
breeder’s desired gain vector **d**. The algorithm in
[`run_desired_gain_index_Joukhadar2024()`](https://FAkohoue.github.io/DGQGSI/reference/run_desired_gain_index_Joukhadar2024.md):

1.  Samples candidate desired-gain vectors in the neighbourhood of
    **d**.
2.  Computes index coefficients **b** for each candidate.
3.  Evaluates the realized gains **g** on the selected subset.
4.  Retains the coefficients that minimise the deviation between **g**
    and **d**.

Because the search is stochastic, multiple independent replicates
(`n_rep`) are run and the best solution across replicates is returned.

### QGSI

For genotype $i$ with GEBV vector ${\mathbf{γ}}_{i}$, the quadratic
index is:

$$\text{QGSI\_DG}_{i} = \mathbf{d}^{\top}{\mathbf{γ}}_{i} + {\mathbf{γ}}_{i}^{\top}\mathbf{W}_{d}\,{\mathbf{γ}}_{i}$$

where **d** is the desired gain vector and $\mathbf{W}_{d}$ is the
quadratic interaction matrix built automatically by
[`construct_Wd_from_dg()`](https://FAkohoue.github.io/DGQGSI/reference/construct_Wd_from_dg.md).
The **linear term** rewards GEBVs aligned with the desired direction;
the **quadratic term** additionally rewards genotypes that combine
favourable trait values simultaneously, capturing synergies and
penalising trade-offs invisible to a linear index.

------------------------------------------------------------------------

## Example data

The package ships with two built-in datasets. `example_pheno.csv`
contains observed phenotypic means per genotype; `example_gebv.csv`
contains genomic estimated breeding values (GEBVs) from a genomic
prediction model.

We work with six traits across two contrasting groups:

| Trait | Description                      | Direction        |
|-------|----------------------------------|------------------|
| YLD   | Grain yield                      | higher is better |
| MY    | Market yield                     | higher is better |
| MI    | Marketable index                 | higher is better |
| BL    | Bacterial leaf blight score      | lower is better  |
| NBL   | Neck bacterial leaf blight score | lower is better  |
| VHB   | Vascular hypersensitive blight   | lower is better  |

The **desired gain vector** `dg` encodes the target improvement per
selection cycle. For the three disease traits, the value represents the
intended *reduction*; direction flipping is managed internally by the
`lower_is_better` argument, so all values in `dg` are given as positive
quantities regardless of the trait direction.

``` r
library(DGQGSI)
library(data.table)

ext   <- system.file("extdata", package = "DGQGSI")
pheno <- fread(file.path(ext, "example_pheno.csv"))
gebv  <- fread(file.path(ext, "example_gebv.csv"))

traits <- c("YLD", "MY", "MI", "BL", "NBL", "VHB")

dg <- c(YLD = 1.5, MY = 0.5, MI = 0.5, BL = 1, NBL = 1, VHB = 1)
```

------------------------------------------------------------------------

## Run the full pipeline

[`run_dgsi_qgsi_pipeline()`](https://FAkohoue.github.io/DGQGSI/reference/run_dgsi_qgsi_pipeline.md)
is the main entry point. Setting `mode = "both"` runs DGSI and QGSI in
sequence and, with `merge_outputs = TRUE`, combines their genotype
rankings into a single comparison table via
[`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DGQGSI/reference/compare_dg_and_qgsi.md).

Key arguments:

| Argument             | Description                                                |
|----------------------|------------------------------------------------------------|
| `mode`               | `"dg"`, `"qgsi"`, or `"both"`                              |
| `init_data`          | Genotype IDs and metadata to preserve in output            |
| `cand_data`          | Candidate phenotypic trait values (DGSI input)             |
| `ref_data`           | Reference population for scaling and covariance estimation |
| `gebv_data`          | Genotype-level GEBVs (QGSI input)                          |
| `trait_cols`         | Character vector of trait column names                     |
| `dg`                 | Named numeric vector of desired gains                      |
| `lower_is_better`    | Traits for which a smaller value is the breeding goal      |
| `trait_min_sd`       | Per-trait minimum threshold for threshold-based selection  |
| `n_iter`             | Optimisation iterations per DGSI replicate                 |
| `n_rep`              | Number of independent DGSI replicates                      |
| `dg_scale_traits`    | Whether to center and scale traits before DGSI             |
| `qgsi_center_traits` | Whether to center GEBVs before QGSI                        |
| `qgsi_scale_traits`  | Whether to scale GEBVs before QGSI                         |
| `merge_outputs`      | Whether to merge DGSI and QGSI rankings into one table     |

``` r
res <- run_dgsi_qgsi_pipeline(
  mode = "both",
  init_data  = pheno[, .(GenoID, Family)],
  cand_data  = pheno[, c("GenoID", traits), with = FALSE],
  ref_data   = pheno[, c("GenoID", traits), with = FALSE],
  gebv_data  = gebv,
  trait_cols = traits,
  dg         = dg,
  lower_is_better = c("BL", "NBL", "VHB"),
  trait_min_sd    = c(YLD = 0.2, MY = 0.1, MI = 0.1,
                      BL  = 0.1, NBL = 0.1, VHB = 0.1),
  dg_scale_traits    = FALSE,
  qgsi_center_traits = FALSE,
  qgsi_scale_traits  = FALSE,
  n_iter        = 50,
  n_rep         = 3,
  merge_outputs = TRUE,
  debug         = FALSE
)
```

The pipeline returns a named list. The components present depend on
`mode`:

- `res$dg_result` — full DGSI output from
  [`run_desired_gain_index_Joukhadar2024()`](https://FAkohoue.github.io/DGQGSI/reference/run_desired_gain_index_Joukhadar2024.md).
- `res$W_d_result` — the automatically constructed $\mathbf{W}_{d}$
  matrix from
  [`construct_Wd_from_dg()`](https://FAkohoue.github.io/DGQGSI/reference/construct_Wd_from_dg.md).
- `res$qgsi_result` — full QGSI output from
  [`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DGQGSI/reference/run_qgsi_desired_gain.md).
- `res$comparison_result` — merged comparison table from
  [`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DGQGSI/reference/compare_dg_and_qgsi.md)
  (only when `mode = "both"` and `merge_outputs = TRUE`).

------------------------------------------------------------------------

## Exploring the results

### The comparison table

The merged comparison table is the primary output for side-by-side
evaluation. It contains, per genotype: the DGSI score
(`DG_SelectionIndex`), the DGSI selection flag (`DG_Selected`), the
linear and quadratic QGSI components (`LinearDGPart`,
`QuadraticDGPart`), the final QGSI score (`QGSI_DG`), the QGSI rank
(`Rank_QGSI_DG`), and a sign-agreement flag (`SignSame_DG_QGSI`).

``` r
comp <- res$comparison_result$comparison_table
head(comp)
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

### Top candidates by DGSI

Sorting by `DG_SelectionIndex` (descending) reveals the genotypes whose
multi-trait phenotypic profile best matches the optimised linear index
weights:

``` r
head(comp[order(-DG_SelectionIndex)], 10)
#>     GenoID DG_SelectionIndex DG_Selected LinearDGPart QuadraticDGPart  QGSI_DG
#>     <char>             <num>      <lgcl>        <num>           <num>    <num>
#>  1:   G013          9.095305        TRUE     3.043771       0.8199342 3.863705
#>  2:   G016          7.375457        TRUE     2.934136       0.7698778 3.704014
#>  3:   G029          7.332733        TRUE     3.245327       0.9101123 4.155439
#>  4:   G001          6.995980        TRUE     2.581310       0.4790381 3.060348
#>  5:   G007          6.072997        TRUE     1.752643       0.2413784 1.994021
#>  6:   G025          5.851942        TRUE     1.849823       1.1382887 2.988112
#>  7:   G006          4.605918        TRUE     1.322001       0.2457242 1.567725
#>  8:   G015          3.600502        TRUE     1.281385       0.4043269 1.685711
#>  9:   G012          3.575645        TRUE     1.849965       0.2775920 2.127557
#> 10:   G003          3.288723        TRUE     1.208304       0.1715159 1.379820
#>     Rank_QGSI_DG SignSame_DG_QGSI
#>            <num>           <lgcl>
#>  1:            2             TRUE
#>  2:            3             TRUE
#>  3:            1             TRUE
#>  4:            4             TRUE
#>  5:            7             TRUE
#>  6:            5             TRUE
#>  7:            9             TRUE
#>  8:            8             TRUE
#>  9:            6             TRUE
#> 10:           10             TRUE
```

### Top candidates by QGSI

Sorting by `Rank_QGSI_DG` (ascending) highlights genotypes that combine
favourable GEBVs across multiple traits simultaneously. The quadratic
term rewards candidates that are good on several fronts at once:

``` r
head(comp[order(Rank_QGSI_DG)], 10)
#>     GenoID DG_SelectionIndex DG_Selected LinearDGPart QuadraticDGPart  QGSI_DG
#>     <char>             <num>      <lgcl>        <num>           <num>    <num>
#>  1:   G029          7.332733        TRUE     3.245327       0.9101123 4.155439
#>  2:   G013          9.095305        TRUE     3.043771       0.8199342 3.863705
#>  3:   G016          7.375457        TRUE     2.934136       0.7698778 3.704014
#>  4:   G001          6.995980        TRUE     2.581310       0.4790381 3.060348
#>  5:   G025          5.851942        TRUE     1.849823       1.1382887 2.988112
#>  6:   G012          3.575645        TRUE     1.849965       0.2775920 2.127557
#>  7:   G007          6.072997        TRUE     1.752643       0.2413784 1.994021
#>  8:   G015          3.600502        TRUE     1.281385       0.4043269 1.685711
#>  9:   G006          4.605918        TRUE     1.322001       0.2457242 1.567725
#> 10:   G003          3.288723        TRUE     1.208304       0.1715159 1.379820
#>     Rank_QGSI_DG SignSame_DG_QGSI
#>            <num>           <lgcl>
#>  1:            1             TRUE
#>  2:            2             TRUE
#>  3:            3             TRUE
#>  4:            4             TRUE
#>  5:            5             TRUE
#>  6:            6             TRUE
#>  7:            7             TRUE
#>  8:            8             TRUE
#>  9:            9             TRUE
#> 10:           10             TRUE
```

### Rank agreement between DGSI and QGSI

A Spearman correlation between the DGSI score ranking and the QGSI rank
summarises overall agreement across all candidates:

``` r
cor_val <- cor(
  rank(-comp$DG_SelectionIndex),
  comp$Rank_QGSI_DG,
  method = "spearman"
)
cat(sprintf("Spearman rank correlation (DGSI vs QGSI): %.3f\n", cor_val))
#> Spearman rank correlation (DGSI vs QGSI): 0.927
```

A high correlation indicates both methods largely agree on which
genotypes to advance. When the correlation is lower, genotypes that
score high on one method but not the other deserve closer inspection:
they tend to have strong trait interactions that only the quadratic term
captures.

The Pearson and Spearman correlations computed by
[`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DGQGSI/reference/compare_dg_and_qgsi.md)
are also available directly:

``` r
res$comparison_result$correlation_summary
#>                      Comparison  Pearson Spearman     N
#>                          <char>    <num>    <num> <int>
#> 1: DG_SelectionIndex vs QGSI_DG 0.921576 0.927475    30
```

------------------------------------------------------------------------

## Inspecting individual workflow outputs

### DGSI output

The DGSI result contains several diagnostics beyond the ranked genotype
table:

``` r
# Optimised index coefficients
res$dg_result$optimized_b
#>        YLD         MY         MI         BL        NBL        VHB 
#>  2.3422139 -1.0119477  0.5052116  1.6742415  3.2816720  1.4745092

# Realized gains in the transformed trait space
res$dg_result$realized_g
#> YLD  MY  MI  BL NBL VHB 
#>   0   0   0   0   0   0

# Objective value of the best replicate (lower = better match to dg)
res$dg_result$best_q
#> [1] 6

# Mean pairwise correlation among replicate index vectors (stability check)
res$dg_result$mean_replicate_cor
#> [1] 0.7722059

# Trait-level summary: all candidates vs. selected subset
res$dg_result$selection_summary
#>     Trait    Mean_All Mean_Selected    SD_All SD_Selected Realized_Gain_SD
#>    <char>       <num>         <num>     <num>       <num>            <num>
#> 1:    YLD  0.49133893    0.49133893 1.0801658   1.0801658                0
#> 2:     MY  0.15544154    0.15544154 0.8077812   0.8077812                0
#> 3:     MI  0.21967789    0.21967789 0.6253735   0.6253735                0
#> 4:     BL  0.14404814    0.14404814 0.9487245   0.9487245                0
#> 5:    NBL -0.03033278   -0.03033278 0.6348214   0.6348214                0
#> 6:    VHB -0.02779384   -0.02779384 0.9537685   0.9537685                0
```

The `selection_summary` table reports the mean and SD for each trait in
the full candidate set and in the selected subset, together with the
realized gain expressed in SD units. Values above 0.9 for
`mean_replicate_cor` indicate good agreement across replicates.

### QGSI output

The QGSI result separates the linear and quadratic components per
genotype:

``` r
head(res$qgsi_result$ranked_geno)
#>    GenoID Family LinearDGPart QuadraticDGPart  QGSI_DG Rank_QGSI_DG
#>    <char> <char>        <num>           <num>    <num>        <num>
#> 1:   G029     F4     3.245327       0.9101123 4.155439            1
#> 2:   G013     F3     3.043771       0.8199342 3.863705            2
#> 3:   G016     F1     2.934136       0.7698778 3.704014            3
#> 4:   G001     F1     2.581310       0.4790381 3.060348            4
#> 5:   G025     F5     1.849823       1.1382887 2.988112            5
#> 6:   G012     F2     1.849965       0.2775920 2.127557            6
#>    Rank_LinearDGPart Rank_QuadraticDGPart
#>                <num>                <num>
#> 1:                 1                    2
#> 2:                 2                    4
#> 3:                 3                    5
#> 4:                 4                    7
#> 5:                 6                    1
#> 6:                 5                   14
```

A component-level summary (mean, SD, range) is also available:

``` r
res$qgsi_result$component_summary
#>          Component      Mean        SD         Min      Max
#>             <char>     <num>     <num>       <num>    <num>
#> 1:    LinearDGPart 0.5528237 1.4480669 -2.65193493 3.245327
#> 2: QuadraticDGPart 0.3666553 0.2782946  0.03645478 1.138289
#> 3:         QGSI_DG 0.9194790 1.5520520 -1.79488272 4.155439
```

------------------------------------------------------------------------

## Running each method separately

Both workflows can also be called independently when finer control is
needed.

### DGSI only

[`run_desired_gain_index_Joukhadar2024()`](https://FAkohoue.github.io/DGQGSI/reference/run_desired_gain_index_Joukhadar2024.md)
supports two selection strategies via `select_mode`:

- `"top_n"` — selects the `n_select` highest-ranking genotypes.
- `"trait_thresholds"` — selects genotypes meeting minimum per-trait
  thresholds (`trait_min_sd`), with optional fallback to top-N if no
  genotype passes.

``` r
dgsi_res <- run_desired_gain_index_Joukhadar2024(
  init_data       = pheno[, .(GenoID, Family)],
  cand_data       = pheno[, c("GenoID", traits), with = FALSE],
  ref_data        = pheno[, c("GenoID", traits), with = FALSE],
  trait_cols      = traits,
  dg              = dg,
  lower_is_better = c("BL", "NBL", "VHB"),
  select_mode     = "trait_thresholds",
  trait_min_sd    = c(YLD = 0.2, MY = 0.1, MI = 0.1,
                      BL  = 0.1, NBL = 0.1, VHB = 0.1),
  scale_traits    = FALSE,
  n_iter          = 200,
  n_rep           = 10,
  debug           = FALSE
)

head(dgsi_res$ranked_geno)
dgsi_res$selection_summary
```

### QGSI only

[`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DGQGSI/reference/run_qgsi_desired_gain.md)
requires GEBVs and a desired-gain vector. A full $\mathbf{W}_{d}$ matrix
can be supplied via `W_d`; otherwise a diagonal matrix is built from
`quadratic_diag_weights`, which captures only squared trait effects
without cross-trait interactions.

``` r
qgsi_res <- run_qgsi_desired_gain(
  init_data         = pheno[, .(GenoID, Family)],
  gebv_data         = gebv[, c("GenoID", traits), with = FALSE],
  trait_cols        = traits,
  dg                = dg,
  lower_is_better   = c("BL", "NBL", "VHB"),
  center_traits     = FALSE,
  scale_traits      = FALSE,
  return_components = TRUE,
  debug             = FALSE
)

head(qgsi_res$ranked_geno)
qgsi_res$component_summary
```

------------------------------------------------------------------------

## The quadratic interaction matrix W_d

[`construct_Wd_from_dg()`](https://FAkohoue.github.io/DGQGSI/reference/construct_Wd_from_dg.md)
derives $\mathbf{W}_{d}$ from the desired gain vector and the empirical
trait correlation structure in the GEBV data. It supports four
construction strategies via the `method` argument:

| Method            | Description                                                           |
|-------------------|-----------------------------------------------------------------------|
| `"hybrid"`        | Combines diagonal, outer-product, and correlation-weighted components |
| `"diag_only"`     | Diagonal matrix only; no cross-trait interactions                     |
| `"dg_outer"`      | Outer product $\mathbf{d}\mathbf{d}^{\top}$ only                      |
| `"corr_weighted"` | Outer product weighted by the empirical trait correlation matrix      |

``` r
Wobj <- construct_Wd_from_dg(
  gebv_data      = gebv[, c("GenoID", traits), with = FALSE],
  trait_cols     = traits,
  dg             = dg,
  lower_is_better = c("BL", "NBL", "VHB"),
  center_traits  = FALSE,
  scale_traits   = FALSE,
  method         = "hybrid",
  lambda_diag    = 1.0,
  lambda_outer   = 0.5,
  lambda_corr    = 1.0,
  offdiag_shrink = 0.5,
  normalize      = "trace",
  debug          = FALSE
)

round(Wobj$W_d, 4)
#>        YLD     MY     MI     BL    NBL    VHB
#> YLD 0.3451 0.0104 0.0135 0.0304 0.0327 0.0277
#> MY  0.0104 0.0619 0.0068 0.0111 0.0156 0.0053
#> MI  0.0135 0.0068 0.0619 0.0118 0.0147 0.0069
#> BL  0.0304 0.0111 0.0118 0.1770 0.0196 0.0118
#> NBL 0.0327 0.0156 0.0147 0.0196 0.1770 0.0115
#> VHB 0.0277 0.0053 0.0069 0.0118 0.0115 0.1770
```

Positive off-diagonal entries indicate that simultaneous improvement in
two traits is rewarded by the quadratic term. The diagonal entries
encode the per-trait squared contribution, scaled to the magnitude of
each desired gain. The empirical trait correlation matrix used to build
$\mathbf{W}_{d}$ is also returned:

``` r
round(Wobj$cor_matrix, 3)
#>        YLD     MY     MI     BL    NBL    VHB
#> YLD  1.000 -0.109  0.007  0.073  0.115  0.023
#> MY  -0.109  1.000  0.269  0.127  0.382 -0.200
#> MI   0.007  0.269  1.000  0.164  0.331 -0.108
#> BL   0.073  0.127  0.164  1.000  0.054 -0.166
#> NBL  0.115  0.382  0.331  0.054  1.000 -0.176
#> VHB  0.023 -0.200 -0.108 -0.166 -0.176  1.000
```

When `mode = "qgsi"` or `mode = "both"` in the pipeline,
$\mathbf{W}_{d}$ is constructed automatically with the same default
settings unless a pre-built matrix is passed via `W_d`.

------------------------------------------------------------------------

## Merging independent outputs

If you run DGSI and QGSI in separate calls,
[`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DGQGSI/reference/compare_dg_and_qgsi.md)
merges them by genotype ID and computes rank-difference diagnostics and
correlation summaries:

``` r
merged <- compare_dg_and_qgsi(
  dg_result                = dgsi_res,
  qgsi_result              = qgsi_res,
  id_col                   = "GenoID",
  compute_rank_differences = TRUE,
  add_correlation_summary  = TRUE,
  sort_by                  = "DG",
  debug                    = FALSE
)

head(merged$comparison_table)
merged$correlation_summary
```

The `sort_by` argument accepts `"DG_rank"`, `"QGSI_rank"`, `"DG"`,
`"QGSI"`, or `"none"`. Note that `"DG_rank"` requires
`Rank_SelectionIndex` to be present in `dg_result$ranked_geno`; use
`"DG"` to sort by score when running `select_mode = "trait_thresholds"`.

------------------------------------------------------------------------

## Practical guidance

**Trait scaling.** If traits are on very different raw scales, set
`dg_scale_traits = TRUE` for DGSI and `qgsi_scale_traits = TRUE` for
QGSI. After scaling, `dg` values are interpreted in standard-deviation
units. When traits are already comparable, leave both flags as `FALSE`
(the default).

**Calibrating desired gains.** The relative magnitudes of values in `dg`
determine the emphasis placed on each trait. Calibrate targets against
historical rates of genetic gain in your germplasm to ensure they are
biologically achievable in a single selection cycle.

**`lower_is_better` traits.** Include every trait for which a decrease
is the breeding goal. Internally, the function multiplies these columns
by $- 1$ so that all objectives point in the same direction. Always
specify `dg` values as positive quantities regardless of trait
direction.

**Optimisation stability.** The DGSI search is stochastic. Increase
`n_rep` and `n_iter` for more stable results at the cost of computation
time. For exploratory work, `n_iter = 50` and `n_rep = 3` are adequate.
For final selection decisions, `n_iter = 200` and `n_rep = 10` or higher
are recommended. The `mean_replicate_cor` diagnostic provides a quick
stability check: values above 0.9 indicate good agreement across
replicates.

**Choosing between DGSI and QGSI.**

| Situation                               | Recommended method                         |
|-----------------------------------------|--------------------------------------------|
| Initial candidate screening             | DGSI — simpler, interpretable coefficients |
| Traits are known to interact            | QGSI — captures nonlinear synergies        |
| Rank stability check desired            | Run both; inspect Spearman correlation     |
| Genomic selection programme with GEBVs  | QGSI — designed for GEBV input             |
| Phenotype-only data, no GEBVs available | DGSI                                       |

**Debug mode.** Set `debug = TRUE` in any function to print intermediate
matrices, gain vectors, and sorting decisions at each step.

------------------------------------------------------------------------

## Citation

Akohoue F (2026). *DGQGSI: Desired-Gain Selection Index and Quadratic
Genomic Selection Index tools for breeding programs*. R package version
0.1.0. <https://github.com/FAkohoue/DGQGSI>

------------------------------------------------------------------------

## References

Joukhadar R, Li Y, Thistlethwaite R, Forrest KL, Tibbits JF, Trethowan
R, Hayden MJ (2024). Optimising desired gain indices to maximise
selection response. *Frontiers in Plant Science* **15**:1337388.
<https://doi.org/10.3389/fpls.2024.1337388>

Cerón-Rojas JJ, Montesinos-López OA, Montesinos-López A, et al. (2026).
Nonlinear genomic selection index accelerates multi-trait crop
improvement. *Nature Communications* **17**:1991.
<https://doi.org/10.1038/s41467-026-69890-3>
