# DGQGSI <img src="man/figures/logo.png" align="right" height="120" alt="DGQGSI logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/FAkohoue/DGQGSI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FAkohoue/DGQGSI/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/FAkohoue/DGQGSI/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/FAkohoue/DGQGSI/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

`DGQGSI` provides tools for two related breeding workflows:

- a desired-gain selection index optimization workflow inspired by Joukhadar et al. (2024)
- a desired-gain adaptation of the quadratic genomic selection index (QGSI)

It also includes:

- automatic construction of the quadratic desired-gain matrix `W_d`
- a full wrapper pipeline for running one or both workflows
- a comparison helper to merge DG and QGSI outputs by genotype

## Included functions

### Core workflows

- `run_desired_gain_index_Joukhadar2024()`
- `run_qgsi_desired_gain()`
- `run_dgsi_qgsi_pipeline()`

### Supporting tools

- `construct_Wd_from_dg()`
- `compare_dg_and_qgsi()`

## Install from local source

```r
install.packages("devtools")
devtools::document("DGQGSI")
devtools::install("DGQGSI", build_vignettes = TRUE)
```

## Install from GitHub

Replace `FAkohoue` with your GitHub username after pushing the repository.

```r
install.packages("remotes")
remotes::install_github("FAkohoue/DGQGSI", build_vignettes = TRUE)
```

## Minimal example

```r
library(DGQGSI)
library(data.table)

traits <- c("YLD", "MY", "MI", "BL", "NBL", "VHB")
dg <- c(YLD = 1.5, MY = 0.5, MI = 0.5, BL = 1, NBL = 1, VHB = 1)

ext <- system.file("extdata", package = "DGQGSI")
pheno <- fread(file.path(ext, "example_pheno.csv"))
gebv  <- fread(file.path(ext, "example_gebv.csv"))

res <- run_dgsi_qgsi_pipeline(
  mode = "both",
  init_data = pheno[, .(GenoID, Family)],
  cand_data = pheno[, c("GenoID", traits), with = FALSE],
  ref_data = pheno[, c("GenoID", traits), with = FALSE],
  gebv_data = gebv[, c("GenoID", traits), with = FALSE],
  trait_cols = traits,
  dg = dg,
  lower_is_better = c("BL", "NBL", "VHB"),
  trait_min_sd = c(YLD = 0.2, MY = 0.1, MI = 0.1, BL = 0.1, NBL = 0.1, VHB = 0.1),
  dg_scale_traits = FALSE,
  qgsi_center_traits = FALSE,
  qgsi_scale_traits = FALSE,
  merge_outputs = TRUE,
  debug = FALSE
)

head(res$comparison_result$comparison_table)
```