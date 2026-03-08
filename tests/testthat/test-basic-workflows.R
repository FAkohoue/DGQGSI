test_that("construct_Wd_from_dg returns square symmetric matrix", {
  ext <- system.file("extdata", package = "DGQGSI")
  gebv <- data.table::fread(file.path(ext, "example_gebv.csv"))
  traits <- c("YLD","MY","MI","BL","NBL","VHB")
  dg <- c(YLD=1.5, MY=0.5, MI=0.5, BL=1, NBL=1, VHB=1)

  obj <- construct_Wd_from_dg(
    gebv_data = gebv,
    trait_cols = traits,
    dg = dg,
    lower_is_better = c("BL", "NBL", "VHB"),
    center_traits = FALSE,
    scale_traits = FALSE,
    debug = FALSE
  )

  expect_true(is.matrix(obj$W_d))
  expect_equal(nrow(obj$W_d), length(traits))
  expect_equal(ncol(obj$W_d), length(traits))
  expect_equal(obj$W_d, t(obj$W_d))
})

test_that("run_qgsi_desired_gain returns ranked output and components", {
  ext <- system.file("extdata", package = "DGQGSI")
  pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
  gebv  <- data.table::fread(file.path(ext, "example_gebv.csv"))
  traits <- c("YLD","MY","MI","BL","NBL","VHB")
  dg <- c(YLD=1.5, MY=0.5, MI=0.5, BL=1, NBL=1, VHB=1)

  res <- run_qgsi_desired_gain(
    init_data = pheno[, .(GenoID, Family)],
    gebv_data = gebv,
    trait_cols = traits,
    dg = dg,
    lower_is_better = c("BL", "NBL", "VHB"),
    center_traits = FALSE,
    scale_traits = FALSE,
    debug = FALSE
  )

  expect_true("ranked_geno" %in% names(res))
  expect_true(all(c("LinearDGPart", "QuadraticDGPart", "QGSI_DG", "Rank_QGSI_DG") %in% names(res$ranked_geno)))
  expect_equal(nrow(res$ranked_geno), nrow(gebv))
})

test_that("run_desired_gain_index_Joukhadar2024 returns selected output", {
  ext <- system.file("extdata", package = "DGQGSI")
  pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
  traits <- c("YLD","MY","MI","BL","NBL","VHB")
  dg <- c(YLD=1.5, MY=0.5, MI=0.5, BL=1, NBL=1, VHB=1)

  res <- run_desired_gain_index_Joukhadar2024(
    init_data = pheno[, .(GenoID, Family)],
    cand_data = pheno[, c("GenoID", traits), with = FALSE],
    ref_data = pheno[, c("GenoID", traits), with = FALSE],
    trait_cols = traits,
    dg = dg,
    lower_is_better = c("BL", "NBL", "VHB"),
    select_mode = "top_n",
    n_select = 5,
    n_iter = 30,
    n_rep = 2,
    debug = FALSE,
    return_all_reps = FALSE
  )

  expect_true("ranked_geno" %in% names(res))
  expect_true(all(c("SelectionIndex", "Selected") %in% names(res$ranked_geno)))
  expect_equal(sum(res$ranked_geno$Selected), 5)
})

test_that("compare_dg_and_qgsi merges outputs correctly", {
  ext <- system.file("extdata", package = "DGQGSI")
  pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
  gebv  <- data.table::fread(file.path(ext, "example_gebv.csv"))
  traits <- c("YLD","MY","MI","BL","NBL","VHB")
  dg <- c(YLD=1.5, MY=0.5, MI=0.5, BL=1, NBL=1, VHB=1)

  res_dg <- run_desired_gain_index_Joukhadar2024(
    init_data = pheno[, .(GenoID, Family)],
    cand_data = pheno[, c("GenoID", traits), with = FALSE],
    ref_data = pheno[, c("GenoID", traits), with = FALSE],
    trait_cols = traits,
    dg = dg,
    lower_is_better = c("BL", "NBL", "VHB"),
    select_mode = "top_n",
    n_select = 5,
    n_iter = 20,
    n_rep = 2,
    debug = FALSE,
    return_all_reps = FALSE
  )

  res_qgsi <- run_qgsi_desired_gain(
    init_data = pheno[, .(GenoID, Family)],
    gebv_data = gebv,
    trait_cols = traits,
    dg = dg,
    lower_is_better = c("BL", "NBL", "VHB"),
    center_traits = FALSE,
    scale_traits = FALSE,
    debug = FALSE
  )

  cmp <- compare_dg_and_qgsi(
    dg_result = res_dg,
    qgsi_result = res_qgsi,
    id_col = "GenoID",
    debug = FALSE
  )

  expect_true("comparison_table" %in% names(cmp))
  expect_true(all(c("DG_SelectionIndex", "QGSI_DG") %in% names(cmp$comparison_table)))
  expect_equal(nrow(cmp$comparison_table), nrow(pheno))
})

test_that("wrapper runs in both mode and returns comparison_result", {
  ext <- system.file("extdata", package = "DGQGSI")
  pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
  gebv  <- data.table::fread(file.path(ext, "example_gebv.csv"))
  traits <- c("YLD","MY","MI","BL","NBL","VHB")
  dg <- c(YLD=1.5, MY=0.5, MI=0.5, BL=1, NBL=1, VHB=1)

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
    n_iter = 20,
    n_rep = 2,
    merge_outputs = TRUE,
    debug = FALSE
  )

  expect_true(all(c("dg_result", "qgsi_result", "comparison_result") %in% names(res)))
})
