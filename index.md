# DGQGSI

![DGQGSI logo](reference/figures/logo.png)

## Overview

`DGQGSI` provides tools for **multi-trait selection in breeding
programs** using desired-gain–based selection indices.

The package integrates two complementary approaches:

- **Desired-Gain Selection Index Optimization (DGSI)**  
  Inspired by Joukhadar et al. (2024), this workflow estimates index
  weights that produce realized genetic gains close to breeder-defined
  desired gains.

- **Desired-Gain Quadratic Genomic Selection Index (QGSI)**  
  Based on Cerón-Rojas et al. (2026), this nonlinear genomic index
  incorporates both linear and quadratic trait components to model
  interactions among traits.

These methods allow breeders to:

- define explicit **desired genetic gains**
- evaluate **multi-trait genomic candidates**
- incorporate **trait interactions**
- compare alternative **selection indices**

The package also includes:

- automatic construction of the quadratic desired-gain matrix `W_d`
- a wrapper pipeline for running one or both workflows
- a comparison helper to merge DGSI and QGSI outputs by genotype

------------------------------------------------------------------------

# Scientific Motivation

Modern crop breeding programs aim to simultaneously improve multiple
traits such as:

- yield
- stress tolerance
- disease resistance
- grain or product quality

Traditional selection indices typically rely on **economic weights**,
which can be difficult to specify and interpret in breeding contexts.

Breeders often prefer to define **desired genetic gains** instead.

DGQGSI provides tools to:

- optimize selection indices using **desired gains**
- integrate **genomic predictions**
- incorporate **trait interactions**
- compare alternative selection strategies

This enables more flexible and biologically meaningful multi-trait
selection.

------------------------------------------------------------------------

# Methodological Background

## 1. Desired-Gain Selection Index Optimization (DGSI)

The desired-gain selection index seeks weights **b** such that the
index:

I = b’X

produces realized gains **g** that match the breeder’s desired gains
**d**.

The algorithm implemented in
[`run_desired_gain_index_Joukhadar2024()`](https://FAkohoue.github.io/DGQGSI/reference/run_desired_gain_index_Joukhadar2024.md):

1.  Samples candidate desired-gain vectors
2.  Computes index weights
3.  Evaluates realized genetic gains
4.  Minimizes deviation between realized and desired gains

This approach follows the framework proposed by:

Joukhadar et al. (2024).

------------------------------------------------------------------------

## 2. Desired-Gain Quadratic Genomic Selection Index (QGSI)

The QGSI extends classical linear indices by incorporating **nonlinear
interactions among traits**.

The index is defined as:

QGSI_i = dᵀ γ_i + γ_iᵀ W_d γ_i

where:

- **d** = desired gain vector  
- **γ_i** = vector of trait GEBVs for genotype i  
- **W_d** = quadratic interaction matrix

The quadratic component allows the index to capture:

- trait synergies
- trade-offs among traits
- nonlinear selection responses

This method is described in:

Cerón-Rojas et al. (2026).

------------------------------------------------------------------------

# Workflow Diagram

The DGQGSI pipeline integrates phenotype or genomic predictions with
desired gain objectives.

    Phenotypes / GEBVs
            │
            │
            ▼
    Desired Gain Vector (d)
            │
            ├───────────────┐
            │               │
            ▼               ▼
    DGSI Optimization     QGSI Computation
            │               │
            ▼               ▼
    Selection Index       QGSI Score
            │               │
            └───────┬───────┘
                    ▼
           Compare Genotype Rankings

------------------------------------------------------------------------

# Main Functions

## Core workflows

- [`run_desired_gain_index_Joukhadar2024()`](https://FAkohoue.github.io/DGQGSI/reference/run_desired_gain_index_Joukhadar2024.md)
- [`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DGQGSI/reference/run_qgsi_desired_gain.md)
- [`run_dgsi_qgsi_pipeline()`](https://FAkohoue.github.io/DGQGSI/reference/run_dgsi_qgsi_pipeline.md)

## Supporting tools

- [`construct_Wd_from_dg()`](https://FAkohoue.github.io/DGQGSI/reference/construct_Wd_from_dg.md)
- [`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DGQGSI/reference/compare_dg_and_qgsi.md)

------------------------------------------------------------------------

# Installation

## Install from GitHub

``` r
install.packages("remotes")
remotes::install_github("FAkohoue/DGQGSI", build_vignettes = TRUE)
```

## Install from local source

``` r
install.packages("devtools")
devtools::install("DGQGSI", build_vignettes = TRUE)
```

------------------------------------------------------------------------

# Minimal Example

``` r
library(DGQGSI)
library(data.table)

traits <- c("YLD", "MY", "MI", "BL", "NBL", "VHB")

dg <- c(
  YLD = 1.5,
  MY  = 0.5,
  MI  = 0.5,
  BL  = 1,
  NBL = 1,
  VHB = 1
)

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
  trait_min_sd = c(
    YLD = 0.2,
    MY  = 0.1,
    MI  = 0.1,
    BL  = 0.1,
    NBL = 0.1,
    VHB = 0.1
  ),
  dg_scale_traits = FALSE,
  qgsi_center_traits = FALSE,
  qgsi_scale_traits = FALSE,
  merge_outputs = TRUE,
  debug = FALSE
)

head(res$comparison_result$comparison_table)
```

------------------------------------------------------------------------

# Documentation

Full documentation and tutorials are available at:

<https://fakohoue.github.io/DGQGSI>

------------------------------------------------------------------------

# Citation

If you use **DGQGSI** in research, please cite the package:

    Akohoue, F. (2026).
    DGQGSI: Desired-Gain Selection Index and Quadratic Genomic Selection Index tools for breeding programs.
    R package version 0.1.0.
    https://github.com/FAkohoue/DGQGSI

------------------------------------------------------------------------

# References

Joukhadar R, Li Y, Thistlethwaite R, Forrest KL, Tibbits JF, Trethowan
R, Hayden MJ (2024).  
Optimising desired gain indices to maximise selection response.  
Frontiers in Plant Science 15:1337388.  
<https://doi.org/10.3389/fpls.2024.1337388>

Cerón-Rojas JJ, Montesinos-López OA, Montesinos-López A, et
al. (2026).  
Nonlinear genomic selection index accelerates multi-trait crop
improvement.  
Nature Communications 17:1991.  
<https://doi.org/10.1038/s41467-026-69890-3>

------------------------------------------------------------------------

# License

MIT License © Félicien Akohoue

------------------------------------------------------------------------

# Author

**Félicien Akohoue**

Alliance Bioversity International & CIAT

Research interests:

- genomic selection
- multi-trait breeding
- quantitative genetics
- crop improvement
