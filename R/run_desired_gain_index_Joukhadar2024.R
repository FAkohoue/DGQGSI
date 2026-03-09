#' Run a desired-gain selection index with optional trait scaling and threshold-based or top-N selection
#'
#' @description
#' Implements a desired-gain selection index optimization strategy inspired by
#' Joukhadar et al. (2024). The function searches for index coefficients that
#' produce realized multi-trait gains close to a user-defined desired-gain
#' vector and then ranks all genotypes using the best replicate.
#'
#' @details
#' The workflow:
#' \enumerate{
#'   \item aligns `init_data` and `cand_data` by `id_col`;
#'   \item optionally flips traits in `lower_is_better` so all objectives are in the same direction;
#'   \item optionally scales traits using `ref_data`;
#'   \item computes a phenotypic correlation matrix `P` and a genetic covariance matrix `G`;
#'   \item repeatedly samples candidate desired-gain vectors, computes index coefficients, evaluates realized gains on the selected set, and keeps the best solution according to a goodness-of-fit criterion.
#' }
#'
#' The function is selection-oriented: in addition to a ranking, it returns the
#' selected subset defined either by `top_n` or by trait thresholds.
#'
#' @param init_data A `data.frame` or `data.table` containing identifiers and metadata to preserve.
#' @param cand_data A `data.frame` or `data.table` containing candidate trait values.
#' @param trait_cols Character vector of trait column names used in the index.
#' @param ref_data Optional reference data used for scaling and for estimating `P` and `G`. If `NULL`, `cand_data` is used.
#' @param id_col Character identifier column. Default is `"GenoID"`.
#' @param scale_traits Logical; if `TRUE`, center and scale traits using `ref_data`. Default is `FALSE`.
#' @param lower_is_better Optional character vector of trait names for which lower values are preferred.
#' @param G Optional square genetic covariance matrix. If `NULL`, estimated from `cov(ref_mat)`.
#' @param dg Named numeric desired-gain vector.
#' @param select_mode One of `"top_n"` or `"trait_thresholds"`.
#' @param n_select Integer number of selected genotypes when `select_mode = "top_n"`.
#' @param trait_min_sd Named numeric vector of minimum thresholds used when `select_mode = "trait_thresholds"`.
#' @param fallback_to_top_n Logical; if `TRUE`, falls back to top-N selection when threshold selection returns zero genotypes.
#' @param n_iter Integer iterations per replicate.
#' @param n_rep Integer number of optimization replicates.
#' @param sd_scale Numeric perturbation scale around the current best desired-gain vector.
#' @param seed Integer random seed.
#' @param ridge_P Numeric ridge penalty added to `P`.
#' @param ridge_M Numeric ridge penalty added to the inner matrix `M`.
#' @param debug Logical; if `TRUE`, print debug messages.
#' @param return_all_reps Logical; if `TRUE`, return all replicate details.
#'
#' @return A list with optimized coefficients, realized gains, replicate diagnostics, ranked genotypes, and selected/non-selected subsets.
#'
#' @import data.table
#' @export
run_desired_gain_index_Joukhadar2024 <- function(
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
    sd_scale = 1.0,
    seed = 42,
    ridge_P = 1e-6,
    ridge_M = 1e-6,
    debug = TRUE,
    return_all_reps = TRUE
) {
  select_mode <- match.arg(select_mode)
  
  .dgqgsi_dbg(debug, "============================================================")
  .dgqgsi_dbg(debug, "Starting run_desired_gain_index_Joukhadar2024()")
  .dgqgsi_dbg(debug, "Selection mode: %s", select_mode)
  .dgqgsi_dbg(debug, "Number of traits: %d", length(trait_cols))
  .dgqgsi_dbg(debug, "Number of replicates: %d | Iterations per replicate: %d", n_rep, n_iter)
  .dgqgsi_dbg(debug, "scale_traits = %s | id_col = %s", scale_traits, id_col)
  
  dt_init <- data.table::as.data.table(data.table::copy(init_data))
  dt_cand <- data.table::as.data.table(data.table::copy(cand_data))
  
  if (is.null(ref_data)) {
    .dgqgsi_dbg(debug, "ref_data is NULL -> using cand_data as reference")
    ref_data <- cand_data
  }
  dt_ref <- data.table::as.data.table(data.table::copy(ref_data))
  
  if (!id_col %in% names(dt_init)) stop(sprintf("id_col '%s' not found in init_data.", id_col), call. = FALSE)
  if (!id_col %in% names(dt_cand)) stop(sprintf("id_col '%s' not found in cand_data.", id_col), call. = FALSE)
  if (!all(trait_cols %in% names(dt_cand))) {
    miss <- setdiff(trait_cols, names(dt_cand))
    stop("Missing trait columns in cand_data: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!all(trait_cols %in% names(dt_ref))) {
    miss <- setdiff(trait_cols, names(dt_ref))
    stop("Missing trait columns in ref_data: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  
  common_ids <- intersect(dt_init[[id_col]], dt_cand[[id_col]])
  if (length(common_ids) == 0L) stop("No common IDs found between init_data and cand_data.", call. = FALSE)
  
  dt_init <- dt_init[get(id_col) %in% common_ids]
  dt_cand <- dt_cand[get(id_col) %in% common_ids]
  data.table::setkeyv(dt_init, id_col)
  data.table::setkeyv(dt_cand, id_col)
  dt_init <- dt_init[dt_cand[[id_col]]]
  
  if (!identical(dt_init[[id_col]], dt_cand[[id_col]])) {
    stop("Row alignment failed: init_data and cand_data IDs are not in identical order after alignment.", call. = FALSE)
  }
  
  for (tr in trait_cols) {
    dt_cand[, (tr) := as.numeric(get(tr))]
    dt_ref[,  (tr) := as.numeric(get(tr))]
  }
  
  for (tr in trait_cols) {
    if (anyNA(dt_cand[[tr]])) {
      mu_tr <- mean(dt_cand[[tr]], na.rm = TRUE)
      dt_cand[is.na(get(tr)), (tr) := mu_tr]
    }
    if (anyNA(dt_ref[[tr]])) {
      mu_tr <- mean(dt_ref[[tr]], na.rm = TRUE)
      dt_ref[is.na(get(tr)), (tr) := mu_tr]
    }
  }
  
  if (!is.null(lower_is_better)) {
    if (!all(lower_is_better %in% trait_cols)) {
      bad <- setdiff(lower_is_better, trait_cols)
      stop("These lower_is_better traits are not in trait_cols: ", paste(bad, collapse = ", "), call. = FALSE)
    }
    dt_cand[, (lower_is_better) := lapply(.SD, function(x) -x), .SDcols = lower_is_better]
    dt_ref[,  (lower_is_better) := lapply(.SD, function(x) -x), .SDcols = lower_is_better]
  }
  
  if (isTRUE(scale_traits)) {
    ref_mat0 <- as.matrix(dt_ref[, trait_cols, with = FALSE])
    mu_ref <- colMeans(ref_mat0)
    sd_ref <- apply(ref_mat0, 2, stats::sd)
    bad_sd <- !is.finite(sd_ref) | sd_ref == 0
    sd_ref[bad_sd] <- 1
    
    scale_with_ref <- function(dt, trait_cols, mu_ref, sd_ref) {
      m <- as.matrix(dt[, trait_cols, with = FALSE])
      m <- sweep(m, 2, mu_ref, FUN = "-")
      m <- sweep(m, 2, sd_ref, FUN = "/")
      out <- data.table::copy(dt)
      out[, (trait_cols) := data.table::as.data.table(m)]
      out
    }
    
    dt_cand <- scale_with_ref(dt_cand, trait_cols, mu_ref, sd_ref)
    dt_ref  <- scale_with_ref(dt_ref, trait_cols, mu_ref, sd_ref)
  }
  
  cand_mat <- as.matrix(dt_cand[, trait_cols, with = FALSE])
  ref_mat  <- as.matrix(dt_ref[, trait_cols, with = FALSE])
  p <- length(trait_cols)
  
  dg <- dg[trait_cols]
  if (anyNA(dg)) stop("dg missing for: ", paste(trait_cols[is.na(dg)], collapse = ", "), call. = FALSE)
  dg <- as.numeric(dg)
  names(dg) <- trait_cols
  
  if (select_mode == "trait_thresholds") {
    if (is.null(trait_min_sd)) stop("trait_min_sd must be provided when select_mode = 'trait_thresholds'.", call. = FALSE)
    trait_min_sd <- trait_min_sd[trait_cols]
    if (anyNA(trait_min_sd)) stop("trait_min_sd missing for: ", paste(trait_cols[is.na(trait_min_sd)], collapse = ", "), call. = FALSE)
    trait_min_sd <- as.numeric(trait_min_sd)
    names(trait_min_sd) <- trait_cols
  }
  
  if (is.null(G)) G <- stats::cov(ref_mat)
  .validate_square_matrix(G, p, "G")
  colnames(G) <- rownames(G) <- trait_cols
  
  P <- stats::cor(ref_mat, use = "pairwise.complete.obs")
  P <- 0.5 * (P + t(P))
  G <- 0.5 * (G + t(G))
  P_r <- P + diag(ridge_P, p)
  
  compute_b <- function(d) {
    X <- solve(P_r, G)
    M <- t(G) %*% X
    M <- 0.5 * (M + t(M)) + diag(ridge_M, p)
    y <- solve(M, matrix(d, ncol = 1))
    as.numeric(X %*% y)
  }
  
  select_rows <- function(idx_vec) {
    if (select_mode == "top_n") {
      ord <- order(idx_vec, decreasing = TRUE)
      return(ord[seq_len(min(n_select, length(ord)))])
    }
    ok <- rep(TRUE, nrow(cand_mat))
    for (j in seq_along(trait_cols)) ok <- ok & (cand_mat[, j] >= trait_min_sd[j])
    keep <- which(ok)
    if (length(keep) == 0L) {
      if (!fallback_to_top_n) stop("No genotypes pass trait_min_sd thresholds.", call. = FALSE)
      ord <- order(idx_vec, decreasing = TRUE)
      keep <- ord[seq_len(min(n_select, length(ord)))]
    }
    keep
  }
  
  realized_g <- function(b) {
    idx <- as.numeric(cand_mat %*% b)
    keep <- select_rows(idx)
    sel <- cand_mat[keep, , drop = FALSE]
    mu_all <- colMeans(cand_mat)
    mu_sel <- colMeans(sel)
    sd_all <- apply(cand_mat, 2, stats::sd)
    sd_all[sd_all == 0 | !is.finite(sd_all)] <- 1
    rg <- (mu_sel - mu_all) / sd_all
    names(rg) <- trait_cols
    rg
  }
  
  gof_one <- function(g_i, dg_i) {
    eps <- 1e-12
    if (sign(g_i) != sign(dg_i)) return(0)
    base <- sqrt(abs(exp(1) * dg_i))
    r <- (base^2) / max(abs(g_i), eps)
    log(r, base = base) / r
  }
  
  q_value <- function(g_vec) {
    gof_i <- mapply(gof_one, g_vec, dg)
    gofMAX_i <- mapply(gof_one, dg, dg)
    gofMAX_i[gofMAX_i == 0] <- 1e-12
    sum((gofMAX_i - gof_i) / gofMAX_i)
  }
  
  summarize_selection <- function(cand_mat, keep, trait_cols) {
    sel_mat <- cand_mat[keep, , drop = FALSE]
    
    mean_all <- colMeans(cand_mat)
    mean_sel <- colMeans(sel_mat)
    sd_all <- apply(cand_mat, 2, stats::sd)
    sd_sel <- apply(sel_mat, 2, stats::sd)
    
    realized <- (mean_sel - mean_all) / sd_all
    realized[!is.finite(sd_all) | sd_all == 0] <- NA_real_
    
    data.table::data.table(
      Trait = trait_cols,
      Mean_All = mean_all,
      Mean_Selected = mean_sel,
      SD_All = sd_all,
      SD_Selected = sd_sel,
      Realized_Gain_SD = realized
    )
  }
  
  set.seed(seed)
  reps <- vector("list", n_rep)
  for (rr in seq_len(n_rep)) {
    mu <- dg
    best_q <- Inf
    best_d <- mu
    best_b <- rep(NA_real_, p)
    best_g <- rep(NA_real_, p)
    q_trace <- numeric(n_iter)
    
    for (it in seq_len(n_iter)) {
      d_samp <- vapply(seq_len(p), function(i) {
        sd_i <- sd_scale * abs(mu[i])
        if (!is.finite(sd_i) || sd_i <= 0) sd_i <- 1e-8
        stats::rnorm(1, mean = mu[i], sd = sd_i)
      }, numeric(1))
      
      b <- compute_b(d_samp)
      g <- realized_g(b)
      q <- q_value(g)
      q_trace[it] <- q
      if (is.finite(q) && q < best_q) {
        best_q <- q
        best_d <- d_samp
        best_b <- b
        best_g <- g
        mu <- d_samp
      }
    }
    
    names(best_d) <- trait_cols
    names(best_b) <- trait_cols
    names(best_g) <- trait_cols
    final_index <- as.numeric(cand_mat %*% best_b)
    keep <- select_rows(final_index)
    ord <- order(final_index, decreasing = TRUE)
    
    reps[[rr]] <- list(
      replicate = rr,
      best_q = best_q,
      best_d = best_d,
      best_b = best_b,
      best_g = best_g,
      index = final_index,
      rank = ord,
      keep = keep,
      q_trace = q_trace
    )
  }
  
  all_index <- do.call(cbind, lapply(reps, `[[`, "index"))
  mean_replicate_cor <- if (ncol(all_index) > 1L) {
    cc <- stats::cor(all_index, use = "pairwise.complete.obs")
    mean(cc[upper.tri(cc)])
  } else NA_real_
  
  best_rep <- which.min(vapply(reps, `[[`, numeric(1), "best_q"))
  best <- reps[[best_rep]]
  
  dt_out <- data.table::copy(dt_init)
  dt_out[, SelectionIndex := best$index]
  dt_out[, Selected := FALSE]
  dt_out[best$keep, Selected := TRUE]
  data.table::setorder(dt_out, -SelectionIndex)
  
  selected_geno <- dt_out[Selected == TRUE]
  non_selected_geno <- dt_out[Selected == FALSE]
  selection_summary <- summarize_selection(cand_mat, best$keep, trait_cols)
  
  out <- list(
    dg = dg,
    trait_min_sd = if (select_mode == "trait_thresholds") trait_min_sd else NULL,
    optimized_d = best$best_d,
    optimized_b = best$best_b,
    realized_g = best$best_g,
    best_q = best$best_q,
    q_trace = best$q_trace,
    mean_replicate_cor = mean_replicate_cor,
    best_replicate = best_rep,
    ranked_geno = dt_out,
    selection_summary = selection_summary,
    selected_geno = selected_geno,
    non_selected_geno = non_selected_geno
  )
  if (isTRUE(return_all_reps)) out$all_reps <- reps
  out
}