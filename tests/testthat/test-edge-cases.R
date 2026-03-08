test_that("construct_Wd_from_dg handles diagonal-only method", {
  ext <- system.file("extdata", package = "DGQGSI")
  gebv <- data.table::fread(file.path(ext, "example_gebv.csv"))
  traits <- c("YLD","MY","MI","BL","NBL","VHB")
  dg <- c(YLD=1.5, MY=0.5, MI=0.5, BL=1, NBL=1, VHB=1)

  obj <- construct_Wd_from_dg(
    gebv_data = gebv,
    trait_cols = traits,
    dg = dg,
    method = "diag_only",
    center_traits = FALSE,
    scale_traits = FALSE,
    debug = FALSE
  )

  expect_equal(obj$W_d[row(obj$W_d) != col(obj$W_d)], rep(0, length(obj$W_d[row(obj$W_d) != col(obj$W_d)])))
})

test_that("run_qgsi_desired_gain accepts user-supplied W_d", {
  ext <- system.file("extdata", package = "DGQGSI")
  pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
  gebv  <- data.table::fread(file.path(ext, "example_gebv.csv"))
  traits <- c("YLD","MY","MI","BL","NBL","VHB")
  dg <- c(YLD=1.5, MY=0.5, MI=0.5, BL=1, NBL=1, VHB=1)
  W_d <- diag(rep(1, length(traits)))
  colnames(W_d) <- rownames(W_d) <- traits

  res <- run_qgsi_desired_gain(
    init_data = pheno[, .(GenoID, Family)],
    gebv_data = gebv,
    trait_cols = traits,
    dg = dg,
    W_d = W_d,
    center_traits = FALSE,
    scale_traits = FALSE,
    debug = FALSE
  )

  expect_true("QGSI_DG" %in% names(res$ranked_geno))
})

test_that("desired-gain function errors when thresholds are missing in threshold mode", {
  ext <- system.file("extdata", package = "DGQGSI")
  pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
  traits <- c("YLD","MY","MI","BL","NBL","VHB")
  dg <- c(YLD=1.5, MY=0.5, MI=0.5, BL=1, NBL=1, VHB=1)

  expect_error(
    run_desired_gain_index_Joukhadar2024(
      init_data = pheno[, .(GenoID, Family)],
      cand_data = pheno[, c("GenoID", traits), with = FALSE],
      ref_data = pheno[, c("GenoID", traits), with = FALSE],
      trait_cols = traits,
      dg = dg,
      select_mode = "trait_thresholds",
      debug = FALSE,
      return_all_reps = FALSE
    ),
    "trait_min_sd must be provided"
  )
})
