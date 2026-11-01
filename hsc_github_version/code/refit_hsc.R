# ============================================================================
# HSC two-branch Tobit dynamic analysis and plotting
# Inputs: confirmed TIC-log1p intensities, fixed trajectory pseudotime, and current metabolite annotations.
# Run: Rscript refit_hsc.R
# Replot only: Rscript refit_hsc.R --plot-only
# ============================================================================

# ===================== 1. Central configuration: usually edit only this section ========================
# Locate the package from this script. When sourcing in RStudio, first set the working directory to the package directory.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
PACKAGE_DIR <- if (length(script_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
} else getwd()
INPUT_DIR <- Sys.getenv("HSC_INPUT_DIR", file.path(PACKAGE_DIR, "input"))

CFG <- list(
  input_intensity = file.path(INPUT_DIR, "log_intensity.csv.gz"),
  input_cells = file.path(INPUT_DIR, "cells.csv"),
  input_annotation = file.path(INPUT_DIR, "annotations.csv"),
  output_dir = Sys.getenv("FIGURE2_HSC_OUTPUT", file.path(PACKAGE_DIR, "runs", format(Sys.time(), "%Y%m%d_%H%M%S_hsc"))),
  # full: refit, test, cluster, and plot from inputs; plot-only: replot archived results.
  run_mode = if ("--plot-only" %in% commandArgs(TRUE)) "plot-only" else "full",
  branches = c("CLP", "LK"),       # Correspond to on_CLP_path/on_LK_path in cells.csv
  censor_limit = 0,                 # Left-censoring limit on the log1p scale, not a measured MS LOD
  spline_df = 3,                    # Natural cubic spline degrees of freedom, tested jointly
  relative_tolerance = 1e-9,       # Record survreg's default relative convergence tolerance explicitly
  cholesky_tolerance = 1e-10,       # Default survreg matrix-decomposition tolerance
  max_iterations = 200,            # Maximum Tobit iterations
  likelihood_tolerance = 1e-5,     # Allowed full/null model numerical error, not a filtering threshold
  prediction_points = 100,         # Curve sample points, not 100 independent time observations
  cores = 4,                       # Use one core on Windows because mclapply is incompatible
  feature_chunk_size = 100,        # Features per batch; affects progress output, not results
  p_adjust_method = "BH",           # Correct all testable feature × branch hypotheses together
  significance_q = 0.05,           # Marks significance; does not filter clustering inputs
  cluster_k = 2:8,                  # Export membership and cluster-mean plotting data for every K
  figure_k = 4,                    # Default plotted K; may be c(4,6) or 2:8
  joint_standardization = "per_branch_then_joint", # Fixed method used by the current report
  # Feature ID highlighted in blue; NULL selects the first available annotated feature.
  selected_feature_id = NULL,
  selected_analysis = "CLP",        # CLP, LK, or joint
  selected_k = 6,
  figure_formats = c("pdf", "png"),# Native R vector PDF and PNG
  png_dpi = 180,
  font_family = "sans",              # Cairo devices support Greek letters and other Unicode annotations
  heatmap_min_limit = 2,
  heatmap_limit_quantile = 0.98,
  profile_width = 13,
  profile_height = 8,
  heatmap_width = 16,
  heatmap_height = 25,
  feature_name_max_chars = 75,
  gray_color = "#77777722",
  mean_color = "#d32f2f",
  selected_color = "#2166ac",
  significant_color = "#b2182b",
  heatmap_colors = c("#2166ac", "#f7f7f7", "#b2182b"),
  seed = Sys.getenv("HSC_RANDOM_SEED", unset = "") # Optional seed for stochastic extensions
)
# Fixed analysis conventions: inputs are already log1p, so do not log again; retain batch; omit cell_type.
# No additional detection-rate, curve-amplitude, or effect-size filters. Invalid models go to audit only.
# ============================================================================

# 2. Dependencies, paths, and input validation --------------------------------
required <- c("survival", "splines", "cluster", "parallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these R packages first: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(cluster))
if (nzchar(CFG$seed)) {
  seed_value <- suppressWarnings(as.integer(CFG$seed))
  if (is.na(seed_value)) stop("HSC_RANDOM_SEED must be a valid integer")
  set.seed(seed_value)
}
if (.Platform$OS.type == "windows") CFG$cores <- 1
stopifnot(CFG$joint_standardization == "per_branch_then_joint",
          CFG$spline_df > 0, CFG$prediction_points > 1,
          all(CFG$figure_k %in% CFG$cluster_k), CFG$selected_k %in% CFG$cluster_k)
OUT <- CFG$output_dir
for (d in c("audit", "model_results", "clusters", "plot_data", "figures")) {
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
}
# CSV retains original feature_id and UTF-8 encoding; prevent R from changing numeric column names to X0.01, etc.
read_table <- function(path) read.csv(path, check.names = FALSE,
                                      stringsAsFactors = FALSE, fileEncoding = "UTF-8")
write_table <- function(x, path) write.csv(x, path, row.names = FALSE,
                                         fileEncoding = "UTF-8", na = "")
cells <- read_table(CFG$input_cells)
annotation <- read_table(CFG$input_annotation)
stopifnot(!anyDuplicated(cells$cell_id), !anyDuplicated(annotation$feature_id))
as_flag <- function(x) tolower(as.character(x)) %in% c("true", "1")
annotation$metabolite <- trimws(gsub("[[:space:]]+", " ", annotation$metabolite))
annotated_ids <- annotation$feature_id[!is.na(annotation$metabolite) &
                                       nzchar(annotation$metabolite)]
write_table(data.frame(input = c(CFG$input_intensity, CFG$input_cells,
                                CFG$input_annotation),
                       md5 = unname(tools::md5sum(c(CFG$input_intensity,
                                 CFG$input_cells, CFG$input_annotation)))),
            file.path(OUT, "audit", "input_md5.csv"))
dput(CFG, file.path(OUT, "audit", paste0("parameters_", CFG$run_mode, ".R")))
capture.output(sessionInfo(), file = file.path(OUT, "audit", "sessionInfo.txt"))

# 3. Tobit models: full and null models ----------------------------------------
fit_one_feature <- function(y, t, batch, grid, weights, id) {
  warnings <- character()
  result <- list(feature_id = id, fit_valid = FALSE, p_value = NA_real_,
    reason = "", n_detected = sum(y > CFG$censor_limit),
    logLik_full = NA_real_, logLik_null = NA_real_, scale = NA_real_,
    curve = rep(NA_real_, length(grid)), observed_curve = rep(NA_real_, length(grid)))
  tryCatch({
    d <- data.frame(y = y, t = t, batch = batch)
    # Surv(...,type='left') treats zero values as censored observations; Gaussian errors give a standard Tobit model.
    full_formula <- Surv(y, y > CFG$censor_limit, type = "left") ~
      batch + ns(t, df = CFG$spline_df)
    null_formula <- Surv(y, y > CFG$censor_limit, type = "left") ~ batch
    fit_model <- function(formula) withCallingHandlers(
      survreg(formula, data = d, dist = "gaussian",
              control = survreg.control(maxiter = CFG$max_iterations,
                rel.tolerance = CFG$relative_tolerance, toler.chol = CFG$cholesky_tolerance)),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    full <- fit_model(full_formula)
    null <- fit_model(null_formula)
    # Apply the current report's numerical validity rules; failure is not the same as a nonsignificant result.
    if (any(!is.finite(coef(full))) || any(!is.finite(coef(null))) ||
        !is.finite(full$scale) || full$scale <= 0 ||
        any(grepl("converge|infinite", warnings, ignore.case = TRUE))) {
      stop("Nonfinite model parameters or nonconvergence: ", paste(unique(warnings), collapse = "; "))
    }
    lf <- as.numeric(logLik(full)); ln <- as.numeric(logLik(null))
    if (!is.finite(lf - ln) || lf < ln - CFG$likelihood_tolerance)
      stop("Full and null log-likelihoods violate the nested-model requirement")
    mu <- observed <- rep(0, length(grid))
    for (i in seq_along(weights)) {
      nd <- data.frame(t = grid, batch = factor(rep(levels(batch)[i], length(grid)),
                                               levels = levels(batch)))
      u <- as.vector(predict(full, nd, type = "lp"))
      mu <- mu + weights[i] * u
      # Export E[max(L,Y*)] as well, but use latent mean mu in the main figure to keep the curves distinct.
      a <- (u - CFG$censor_limit) / full$scale
      expected <- CFG$censor_limit * pnorm(-a) + u * pnorm(a) + full$scale * dnorm(a)
      observed <- observed + weights[i] * expected
    }
    if (any(!is.finite(mu))) stop("Predicted curve contains nonfinite values")
    result$fit_valid <- TRUE
    result$p_value <- pchisq(max(0, 2 * (lf - ln)), CFG$spline_df, lower.tail = FALSE)
    result$logLik_full <- lf; result$logLik_null <- ln; result$scale <- full$scale
    result$curve <- mu; result$observed_curve <- observed
    result$reason <- paste(unique(warnings), collapse = "; ")
    result
  }, error = function(e) { result$reason <- conditionMessage(e); result })
}

if (CFG$run_mode == "full") {
  intensity <- read_table(CFG$input_intensity)
  stopifnot(!anyDuplicated(intensity$cell_id),
            setequal(intensity$cell_id, cells$cell_id))
  intensity <- intensity[match(cells$cell_id, intensity$cell_id), , drop = FALSE]
  Y <- as.matrix(intensity[, setdiff(names(intensity), "cell_id"), drop = FALSE])
  storage.mode(Y) <- "double"
  if (any(!is.finite(Y)) || any(Y < CFG$censor_limit))
    stop("Input must contain finite TIC-log1p intensities no lower than the specified censoring limit")
  ids <- colnames(Y); stopifnot(all(ids %in% annotation$feature_id))
  grid <- seq(0, 1, length.out = CFG$prediction_points)
  fits <- stats <- normalized_cells <- list()
  for (br in CFG$branches) {
    ix <- as_flag(cells[[paste0("on_", br, "_path")]])
    t <- cells$Pseudotime[ix]; if (diff(range(t)) <= 0) stop("Pseudotime range is zero")
    t <- (t - min(t)) / diff(range(t))
    batch <- factor(cells$batch[ix]); weights <- as.numeric(table(batch)) / length(batch)
    A <- Y[ix, , drop = FALSE]
    normalized_cells[[br]] <- cbind(branch = br, cells[ix, , drop = FALSE],
                                    normalized_pseudotime = t)
    write_table(data.frame(batch = levels(batch), weight = weights),
                file.path(OUT, "audit", paste0(br, "_prediction_batch_weights.csv")))
    rr <- list()
    for (start in seq(1, ncol(A), by = CFG$feature_chunk_size)) {
      jj <- start:min(ncol(A), start + CFG$feature_chunk_size - 1)
      part <- parallel::mclapply(jj, function(j)
        fit_one_feature(A[, j], t, batch, grid, weights, ids[j]), mc.cores = CFG$cores)
      rr <- c(rr, part)
      message(br, ": completed ", max(jj), "/", ncol(A), " features")
    }
    fits[[br]] <- rr
    stats[[br]] <- do.call(rbind, lapply(rr, function(r)
      as.data.frame(r[!names(r) %in% c("curve", "observed_curve")], stringsAsFactors = FALSE)))
  }
  # Correction covers all testable hypotheses in both branches; do not replace NA with p=0 or p=1.
  all_p <- unlist(lapply(stats, `[[`, "p_value"), use.names = FALSE)
  all_q <- rep(NA_real_, length(all_p)); good <- is.finite(all_p)
  all_q[good] <- p.adjust(all_p[good], method = CFG$p_adjust_method)
  offset <- 0
  for (br in CFG$branches) {
    n <- nrow(stats[[br]])
    stats[[br]]$q_value <- all_q[offset + seq_len(n)]; offset <- offset + n
    stats[[br]]$branch_q <- p.adjust(stats[[br]]$p_value, method = CFG$p_adjust_method)
    st <- stats[[br]]; valid <- st$fit_valid
    write_table(st, file.path(OUT, "audit", paste0(br, "_all_fit_statistics.csv")))
    write_table(st[!valid, ], file.path(OUT, "audit", paste0(br, "_excluded_features.csv")))
    write_table(st[valid, ], file.path(OUT, "model_results", paste0(br, "_statistics.csv")))
    for (kind in c("curve", "observed_curve")) {
      mat <- do.call(rbind, lapply(fits[[br]], `[[`, kind))
      colnames(mat) <- format(grid, digits = 12)
      suffix <- if (kind == "curve") "predicted_log_curves" else "expected_observed_log_curves"
      write_table(data.frame(feature_id = ids[valid], mat[valid, , drop = FALSE], check.names = FALSE),
                  file.path(OUT, "model_results", paste0(br, "_", suffix, ".csv")))
    }
  }
  write_table(do.call(rbind, normalized_cells), file.path(OUT, "plot_data", "cell_pseudotime.csv"))

  # 4. Cluster all valid curves; significance does not filter clustering inputs --------
  mats <- list()
  for (br in CFG$branches) {
    d <- read_table(file.path(OUT, "model_results", paste0(br, "_predicted_log_curves.csv")))
    M <- as.matrix(d[, -1]); rownames(M) <- d$feature_id
    usable <- apply(M, 1, function(v) all(is.finite(v)) && sd(v) > 0)
    mats[[br]] <- M[usable, , drop = FALSE]
  }
  for (br in c(CFG$branches, "joint")) {
    if (br == "joint") {
      ids <- intersect(rownames(mats[[1]]), rownames(mats[[2]]))
      # Standardize branches separately before concatenation to compare shapes, not arbitrary latent baselines.
      X <- cbind(t(scale(t(mats[[1]][ids, ]))), t(scale(t(mats[[2]][ids, ]))))
    } else { X <- mats[[br]]; ids <- rownames(X) }
    Z <- t(scale(t(X)))                 # Use the sample SD, consistent with the current report
    distance <- pmax((1 - cor(t(Z))) / 2, 0); diag(distance) <- 0
    diagnostic <- list()
    for (k in CFG$cluster_k) {
      if (k >= nrow(Z)) stop("Cluster count must be smaller than the number of valid features")
      fit <- cluster::pam(as.dist(distance), k, diss = TRUE)
      member <- data.frame(feature_id = ids, cluster = as.integer(fit$clustering),
                           is_medoid = seq_along(ids) %in% fit$id.med)
      ai <- match(ids, annotation$feature_id)
      member$metabolite <- annotation$metabolite[ai]
      member$metabolite_class <- annotation$metabolite_class[ai]
      member$mz_center <- annotation$mz_center[ai]
      for (b in CFG$branches) {
        st <- stats[[b]]; member[[paste0(b, "_q_value")]] <- st$q_value[match(ids, st$feature_id)]
      }
      write_table(member, file.path(OUT, "clusters", paste0(br, "_K", k, "_members.csv")))
      means <- sapply(seq_len(k), function(j) colMeans(Z[fit$clustering == j, , drop = FALSE]))
      cr <- cor(means); diag(cr) <- NA
      diagnostic[[as.character(k)]] <- data.frame(k = k, silhouette = fit$silinfo$avg.width,
        min_cluster_size = min(table(fit$clustering)), max_cluster_mean_correlation = max(cr, na.rm = TRUE))
    }
    # Column order is documented in prediction_grid.csv; the joint matrix has 100 CLP points, then 100 LK points.
    write_table(data.frame(feature_id = ids, Z, check.names = FALSE),
                file.path(OUT, "clusters", paste0(br, "_standardized_curves.csv")))
    write_table(do.call(rbind, diagnostic), file.path(OUT, "clusters", paste0(br, "_k_diagnostics.csv")))
  }
} else if (CFG$run_mode != "plot-only") stop("run_mode must be full or plot-only")

# 5. Plotting utilities and standalone plotting-data exports ------------------
# Subsequent steps read only valid membership tables and curve matrices; untestable features never enter figures or selection lists.
read_curves <- function(br) {
  d <- read_table(file.path(OUT, "clusters", paste0(br, "_standardized_curves.csv")))
  M <- as.matrix(d[, -1]); rownames(M) <- d$feature_id; M
}
read_members <- function(br, k) read_table(file.path(OUT, "clusters", paste0(br, "_K", k, "_members.csv")))
grid <- seq(0, 1, length.out = CFG$prediction_points)
write_table(data.frame(column_index = seq_len(2 * length(grid)),
                       branch = rep(CFG$branches, each = length(grid)),
                       normalized_pseudotime = rep(grid, 2)),
            file.path(OUT, "plot_data", "joint_prediction_grid.csv"))

# Draw the same content through a callback on PDF and PNG devices to avoid inconsistent output.
save_figure <- function(stem, width, height, draw) {
  for (fmt in CFG$figure_formats) {
    path <- file.path(OUT, "figures", paste0(stem, ".", fmt))
    if (fmt == "pdf") {
      if (capabilities("cairo")) cairo_pdf(path, width = width, height = height, family = CFG$font_family)
      else pdf(path, width = width, height = height, useDingbats = FALSE)
    }
    else if (fmt == "png") png(path, width = width, height = height, units = "in", res = CFG$png_dpi)
    else stop("Only pdf and png formats are currently supported")
    tryCatch(draw(), finally = dev.off())
  }
}
draw_curve_panel <- function(values, title, selected = NULL, selected_label = "Medoid") {
  matplot(grid, t(values), type = "l", lty = 1, col = CFG$gray_color, lwd = 0.5,
          xlab = "Normalized pseudotime", ylab = "Standardized latent log1p intensity", main = title)
  lines(grid, colMeans(values), col = CFG$mean_color, lwd = 2.5)
  if (!is.null(selected)) lines(grid, selected, col = CFG$selected_color, lwd = 2)
  legend("topleft", legend = c("Cluster mean", selected_label),
         col = c(CFG$mean_color, CFG$selected_color), lty = 1, lwd = 2, bty = "n", cex = 0.65)
}

all_means <- list(); entry <- 0
for (br in c(CFG$branches, "joint")) {
  Z <- read_curves(br); paths <- if (br == "joint") CFG$branches else br
  for (k in CFG$cluster_k) {
    member <- read_members(br, k)
    annotated <- member[member$feature_id %in% annotated_ids, , drop = FALSE]
    annotated <- annotated[order(annotated$cluster, annotated$metabolite, annotated$feature_id), ]
    for (b in CFG$branches) annotated[[paste0(b, "_significant")]] <-
      !is.na(annotated[[paste0(b, "_q_value")]]) & annotated[[paste0(b, "_q_value")]] < CFG$significance_q
    write_table(annotated, file.path(OUT, "plot_data", paste0(br, "_K", k, "_annotated.csv")))
    for (cl in sort(unique(member$cluster))) {
      g <- member[member$cluster == cl, ]; V <- Z[g$feature_id, , drop = FALSE]
      for (j in seq_along(paths)) {
        cols <- (j - 1) * length(grid) + seq_along(grid); entry <- entry + 1
        all_means[[entry]] <- data.frame(analysis = br, K = k, cluster = cl, branch = paths[j],
          normalized_pseudotime = grid, mean_standardized = colMeans(V[, cols, drop = FALSE]),
          n_features = nrow(V))
      }
    }
    if (!k %in% CFG$figure_k) next
    save_figure(paste0(br, "_K", k, "_profiles"), CFG$profile_width,
                if (br == "joint") 2.5 * k else CFG$profile_height, function() {
      par(mfrow = if (br == "joint") c(k, 2) else c(ceiling(k / 3), 3),
          mar = c(4, 4, 3, 1), cex = 0.8)
      for (cl in sort(unique(member$cluster))) {
        g <- member[member$cluster == cl, ]; medoid <- g$feature_id[as_flag(g$is_medoid)][1]
        for (j in seq_along(paths)) {
          cols <- (j - 1) * length(grid) + seq_along(grid)
          draw_curve_panel(Z[g$feature_id, cols, drop = FALSE],
            paste(paths[j], "C", cl, " n=", nrow(g)), Z[medoid, cols])
        }
      }
    })
    # The long heatmap table records feature, cluster, pseudotime, standardized value, and star criteria for replotting elsewhere.
    heat_long <- lapply(seq_along(paths), function(j) {
      cols <- (j - 1) * length(grid) + seq_along(grid)
      data.frame(analysis = br, K = k, branch = paths[j],
        feature_id = rep(annotated$feature_id, each = length(grid)),
        metabolite = rep(annotated$metabolite, each = length(grid)),
        cluster = rep(annotated$cluster, each = length(grid)),
        row_order = rep(seq_len(nrow(annotated)), each = length(grid)),
        normalized_pseudotime = rep(grid, nrow(annotated)),
        standardized_value = as.vector(t(Z[annotated$feature_id, cols, drop = FALSE])),
        q_value = rep(annotated[[paste0(paths[j], "_q_value")]], each = length(grid)),
        significant = rep(annotated[[paste0(paths[j], "_significant")]], each = length(grid)))
    })
    write_table(do.call(rbind, heat_long), file.path(OUT, "plot_data", paste0(br, "_K", k, "_annotated_heatmap_long.csv")))
    save_figure(paste0(br, "_K", k, "_annotated_heatmap"),
                CFG$heatmap_width * length(paths), CFG$heatmap_height, function() {
      par(mfrow = c(1, length(paths)), mar = c(4, 24, 4, 4))
      for (j in seq_along(paths)) {
        cols <- (j - 1) * length(grid) + seq_along(grid)
        V <- Z[annotated$feature_id, cols, drop = FALSE]; n <- nrow(V)
        lim <- max(CFG$heatmap_min_limit, as.numeric(quantile(abs(V), CFG$heatmap_limit_quantile)))
        pal <- colorRampPalette(CFG$heatmap_colors)(101)
        Vplot <- V; Vplot[Vplot > lim] <- lim; Vplot[Vplot < -lim] <- -lim
        image(grid, seq_len(n), t(Vplot[n:1, , drop = FALSE]), col = pal,
              zlim = c(-lim, lim), axes = FALSE, xlab = "Normalized pseudotime", ylab = "",
              main = paste(paths[j], "K =", k, "|", n, "annotated features"))
        axis(1, at = seq(0, 1, .25)); box()
        for (i in seq_len(n)) {
          sig <- annotated[[paste0(paths[j], "_significant")]][i]
          lab <- substr(annotated$metabolite[i], 1, CFG$feature_name_max_chars)
          axis(2, at = n - i + 1, labels = paste0(if (sig) "* " else "", lab),
               las = 2, tick = FALSE, cex.axis = .55,
               col.axis = if (sig) CFG$significant_color else "black", font = if (sig) 2 else 1)
        }
        bounds <- cumsum(table(factor(annotated$cluster, levels = sort(unique(annotated$cluster)))))
        abline(h = n - head(bounds, -1) + .5, col = "black", lwd = .7)
        centers <- n - (c(0, head(bounds, -1)) + bounds) / 2 + .5
        axis(4, at = centers, labels = paste0("C", sort(unique(annotated$cluster))),
             las = 2, tick = FALSE, cex.axis = .75)
        mtext(paste0("* BH q < ", CFG$significance_q,
                     "; colors: blue=-", round(lim, 2), ", white=0, red=+", round(lim, 2)),
              side = 3, line = .3, cex = .65)
      }
    })
  }
}
write_table(do.call(rbind, all_means), file.path(OUT, "plot_data", "cluster_mean_curves_long.csv"))

# 6. Selected annotated feature: blue is the feature, red is its cluster-wide mean --------
br <- CFG$selected_analysis; k <- CFG$selected_k
Z <- read_curves(br); member <- read_members(br, k)
available <- member$feature_id[member$feature_id %in% annotated_ids]
fid <- if (is.null(CFG$selected_feature_id)) available[1] else CFG$selected_feature_id
if (!fid %in% available) stop("The selected feature_id has no valid annotated curve in this analysis")
cl <- member$cluster[match(fid, member$feature_id)]
group_ids <- member$feature_id[member$cluster == cl]
paths <- if (br == "joint") CFG$branches else br
name <- annotation$metabolite[match(fid, annotation$feature_id)]
selected_data <- lapply(seq_along(paths), function(j) {
  cols <- (j - 1) * length(grid) + seq_along(grid)
  data.frame(analysis = br, K = k, branch = paths[j], feature_id = fid, metabolite = name,
    cluster = cl, normalized_pseudotime = grid, selected_value = Z[fid, cols],
    cluster_mean = colMeans(Z[group_ids, cols, drop = FALSE]))
})
write_table(do.call(rbind, selected_data), file.path(OUT, "plot_data", "selected_feature_curves.csv"))
save_figure("selected_feature", 8 * length(paths), 6, function() {
  par(mfrow = c(1, length(paths)), mar = c(7, 4, 3, 1))
  for (j in seq_along(paths)) {
    cols <- (j - 1) * length(grid) + seq_along(grid)
    draw_curve_panel(Z[group_ids, cols, drop = FALSE], paste(paths[j], "Cluster", cl),
                     Z[fid, cols], "Selected feature")
    mtext(paste(strwrap(name, width = 70), collapse = "\n"), side = 1, line = 4,
          col = CFG$selected_color, cex = .8)
  }
})

# The cross-branch table retains only features valid on both sides so untested and nonsignificant are not conflated.
left <- read_members(CFG$branches[1], CFG$selected_k)[, c("feature_id", "cluster")]
right <- read_members(CFG$branches[2], CFG$selected_k)[, c("feature_id", "cluster")]
names(left)[2] <- paste0(CFG$branches[1], "_cluster")
names(right)[2] <- paste0(CFG$branches[2], "_cluster")
write_table(merge(left, right, by = "feature_id"),
            file.path(OUT, "plot_data", "branch_cluster_correspondence.csv"))
message("Analysis and plotting complete. Output directory: ", normalizePath(OUT))
