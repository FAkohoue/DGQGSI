#' Run a desired-gain selection index with optional trait scaling and
#' threshold-based or top-N selection
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
#'   \item optionally flips traits in `lower_is_better` so all objectives are in
#'         the same direction;
#'   \item optionally scales traits using `ref_data`;
#'   \item computes a phenotypic correlation matrix `P` and a genetic covariance
#'         matrix `G`;
#'   \item repeatedly samples candidate desired-gain vectors, computes index
#'         coefficients, evaluates realized gains on the selected set, and keeps
#'         the best solution according to a goodness-of-fit criterion.
#' }
#'
#' This function is **selection-oriented**. It does not only return one index
#' value per genotype, but also determines a selected subset based on either:
#' \itemize{
#'   \item the highest index values (`select_mode = "top_n"`), or
#'   \item user-defined minimum trait thresholds
#'         (`select_mode = "trait_thresholds"`).
#' }
#'
#' If your trait values are already expressed in comparable standardized units
#' (for example SD units), set `scale_traits = FALSE`. If raw trait values are
#' used, `scale_traits = TRUE` is generally more appropriate.
#'
#' @param init_data A `data.frame` or `data.table` containing genotype
#'   identifiers and any metadata you want to preserve in the final output.
#'   Typical columns may include `GenoID`, family, group, or other descriptive
#'   information. The function returns this table augmented with the final
#'   `SelectionIndex` and `Selected` columns.
#'
#' @param cand_data A `data.frame` or `data.table` containing the candidate
#'   trait values used to build the desired-gain index. It must contain
#'   `id_col` and all columns listed in `trait_cols`. Each row represents one
#'   genotype and each trait column should be numeric or coercible to numeric.
#'
#' @param trait_cols Character vector giving the names of the trait columns to
#'   include in the index. These names must be present in both `cand_data` and
#'   `ref_data` (or in `cand_data` alone if `ref_data = NULL`). The order of
#'   `trait_cols` defines the order used for `dg`, `trait_min_sd`, `P`, and `G`.
#'
#' @param ref_data Optional `data.frame` or `data.table` used as the reference
#'   population for scaling traits and for estimating the phenotypic correlation
#'   matrix `P` and, if necessary, the genetic covariance matrix `G`. If `NULL`,
#'   `cand_data` is used as the reference.
#'
#' @param id_col Character string naming the genotype identifier column used to
#'   align `init_data` and `cand_data`. Default is `"GenoID"`. The function
#'   requires that this column exists in both `init_data` and `cand_data`.
#'
#' @param scale_traits Logical. If `TRUE`, the trait columns in both candidate
#'   and reference data are centered and scaled using the mean and standard
#'   deviation of `ref_data`. Use this when traits are on different raw scales.
#'   If `FALSE`, the function assumes traits are already on a comparable scale,
#'   such as SD units.
#'
#' @param lower_is_better Optional character vector naming traits for which
#'   smaller values are favorable. These traits are internally multiplied by
#'   `-1` so that all traits are oriented in a common “higher is better”
#'   direction during optimization.
#'
#' @param G Optional square genetic covariance matrix among traits. Its dimension
#'   must be `length(trait_cols) x length(trait_cols)`. If `NULL`, the function
#'   estimates `G` as `cov(ref_mat)` from the reference trait matrix.
#'
#' @param dg Named numeric vector of desired gains, with names matching
#'   `trait_cols`. Each value represents the desired response for one trait in
#'   the trait space used by the function. For best interpretability, `dg`
#'   should be on the same scale as the transformed trait values.
#'
#' @param select_mode Character string specifying how the selected subset is
#'   determined. One of:
#'   \itemize{
#'     \item `"top_n"`: selects the `n_select` highest-ranked genotypes by the
#'           final index;
#'     \item `"trait_thresholds"`: selects genotypes meeting minimum trait
#'           thresholds in the transformed trait space.
#'   }
#'
#' @param n_select Integer giving the number of genotypes to select when
#'   `select_mode = "top_n"`. Ignored when `select_mode = "trait_thresholds"`
#'   unless fallback is triggered.
#'
#' @param trait_min_sd Named numeric vector of minimum acceptable trait values
#'   used when `select_mode = "trait_thresholds"`. Names must match
#'   `trait_cols`. Thresholds are applied after direction flipping and optional
#'   scaling, so they should be specified on that transformed scale.
#'
#' @param fallback_to_top_n Logical. Relevant only when
#'   `select_mode = "trait_thresholds"`. If `TRUE` and no genotype passes the
#'   thresholds, the function falls back to selecting the top `n_select`
#'   genotypes by the final index. If `FALSE`, the function stops with an error.
#'
#' @param n_iter Integer giving the number of optimization iterations performed
#'   within each replicate. Larger values increase the search effort but also
#'   increase computation time.
#'
#' @param n_rep Integer giving the number of independent optimization replicates.
#'   Replicates are useful because the search is stochastic, and different
#'   replicates may converge to slightly different solutions.
#'
#' @param sd_scale Numeric factor controlling the perturbation scale used when
#'   sampling candidate desired-gain vectors around the current best solution.
#'   Larger values increase exploration; smaller values keep the search more
#'   local around the current optimum.
#'
#' @param seed Integer random seed used to make the stochastic optimization
#'   reproducible.
#'
#' @param ridge_P Numeric ridge penalty added to the phenotypic correlation
#'   matrix `P` before inversion. This can improve numerical stability when `P`
#'   is close to singular.
#'
#' @param ridge_M Numeric ridge penalty added to the inner matrix `M` during the
#'   coefficient computation. This can improve stability in cases of highly
#'   correlated traits or poorly conditioned matrices.
#'
#' @param debug Logical. If `TRUE`, the function prints progress and diagnostic
#'   messages during execution. This is helpful for understanding what the
#'   function is doing and for troubleshooting.
#'
#' @param return_all_reps Logical. If `TRUE`, the returned object includes the
#'   full list of replicate-level optimization results. If `FALSE`, only the
#'   summary corresponding to the best replicate is returned.
#'
#' @return A list containing:
#' \describe{
#'   \item{dg}{The desired-gain vector used in the analysis.}
#'   \item{trait_min_sd}{The threshold vector used when
#'     `select_mode = "trait_thresholds"`, otherwise `NULL`.}
#'   \item{optimized_d}{The best sampled desired-gain vector from the winning
#'     replicate.}
#'   \item{optimized_b}{The optimized index coefficient vector from the winning
#'     replicate.}
#'   \item{realized_g}{The realized gains achieved by the selected subset in the
#'     transformed trait space.}
#'   \item{best_q}{The objective value of the best replicate. Lower values
#'     indicate a closer match between realized and desired gains.}
#'   \item{q_trace}{The iteration-wise objective values from the best replicate.}
#'   \item{mean_replicate_cor}{Mean pairwise correlation among replicate index
#'     vectors, used as a simple replicate-stability diagnostic.}
#'   \item{best_replicate}{The index of the replicate with the best objective
#'     value.}
#'   \item{ranked_geno}{A `data.table` containing `init_data` plus the final
#'     `SelectionIndex` and `Selected` columns, sorted by decreasing
#'     `SelectionIndex`.}
#'   \item{selection_summary}{A `data.table` summarizing mean trait values and
#'     realized gains for all candidates versus the selected subset.}
#'   \item{selected_geno}{Subset of `ranked_geno` with `Selected == TRUE`.}
#'   \item{non_selected_geno}{Subset of `ranked_geno` with `Selected == FALSE`.}
#'   \item{all_reps}{Optional list of all replicate results when
#'     `return_all_reps = TRUE`.}
#' }
#'
#' @examples
#' trait_cols <- c("YLD", "MY", "MI", "BL", "NBL", "VHB")
#'
#' dg <- c(
#'   YLD = 1.5,
#'   MY  = 0.5,
#'   MI  = 0.5,
#'   BL  = 1.0,
#'   NBL = 1.0,
#'   VHB = 1.0
#' )
#'
#' trait_min_sd <- c(
#'   YLD = 0.2,
#'   MY  = 0.1,
#'   MI  = 0.1,
#'   BL  = 0.1,
#'   NBL = 0.1,
#'   VHB = 0.1
#' )
#'
#' ext <- system.file("extdata", package = "DGQGSI")
#' pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
#'
#' res <- run_desired_gain_index_Joukhadar2024(
#'   init_data = pheno[, .(GenoID, Family)],
#'   cand_data = pheno[, c("GenoID", trait_cols), with = FALSE],
#'   ref_data = pheno[, c("GenoID", trait_cols), with = FALSE],
#'   trait_cols = trait_cols,
#'   dg = dg,
#'   lower_is_better = c("BL", "NBL", "VHB"),
#'   select_mode = "trait_thresholds",
#'   trait_min_sd = trait_min_sd,
#'   scale_traits = FALSE,
#'   debug = FALSE
#' )
#'
#' head(res$ranked_geno)
#' res$selection_summary
#'
#' @references
#' Joukhadar R, Li Y, Thistlethwaite R, Forrest KL, Tibbits JF, Trethowan R,
#' Hayden MJ (2024). Optimising desired gain indices to maximise selection
#' response. \emph{Frontiers in Plant Science}, 15:1337388.
#' doi:10.3389/fpls.2024.1337388
#'
#' @import data.table
#' @export
#' 
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