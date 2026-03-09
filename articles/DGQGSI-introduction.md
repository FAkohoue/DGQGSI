# DGQGSI introduction

DGQGSI provides two complementary workflows:

- a desired-gain selection index optimization function
- a desired-gain adaptation of QGSI based on trait GEBVs

## Example data

``` r
ext <- system.file("extdata", package = "DGQGSI")
pheno <- fread(file.path(ext, "example_pheno.csv"))
gebv  <- fread(file.path(ext, "example_gebv.csv"))
traits <- c("YLD","MY","MI","BL","NBL","VHB")
dg <- c(YLD=1.5, MY=0.5, MI=0.5, BL=1, NBL=1, VHB=1)
```

## Run the full pipeline

``` r
res <- run_dgsi_qgsi_pipeline(
  mode = "both",
  init_data = pheno[, .(GenoID, Family)],
  cand_data = pheno[, c("GenoID", traits), with = FALSE],
  ref_data = pheno[, c("GenoID", traits), with = FALSE],
  gebv_data = gebv,
  trait_cols = traits,
  dg = dg,
  lower_is_better = c("BL", "NBL", "VHB"),
  trait_min_sd = c(YLD=0.2, MY=0.1, MI=0.1, BL=0.1, NBL=0.1, VHB=0.1),
  dg_scale_traits = FALSE,
  qgsi_center_traits = FALSE,
  qgsi_scale_traits = FALSE,
  n_iter = 50,
  n_rep = 3,
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
