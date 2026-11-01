# MPP explorer: scope-specific Tobit models ----
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(splines))

find_project <- function() {
  candidates <- c(getwd())
  for (f in sys.frames()) if (!is.null(f$ofile)) candidates <- c(candidates, dirname(f$ofile))
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(arg)) candidates <- c(candidates, dirname(sub("^--file=", "", arg[1])))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    candidates <- c(candidates, dirname(rstudioapi::getActiveDocumentContext()$path))
  for (p in candidates) if (file.exists(file.path(p, "input/mpp_heatmap/cells.csv")))
    return(normalizePath(p, winslash = "/"))
  stop("Set FIGURE2_ROOT to the external data-bundle root.")
}
project <- Sys.getenv("FIGURE2_ROOT"); if (!nzchar(project) && !nzchar(Sys.getenv("FIGURE2_INPUT_DIR"))) project <- find_project()
input <- Sys.getenv("FIGURE2_INPUT_DIR", file.path(project, "input/mpp_heatmap"))
output <- Sys.getenv("FIGURE2_MPP_OUTPUT", file.path(project, "runs", format(Sys.time(), "%Y%m%d_%H%M%S_mpp")))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
meta <- read.csv(file.path(input, "cells.csv"))
dat <- read.csv(file.path(input, "annotated_log_intensity.csv"), check.names = FALSE)
stopifnot(identical(dat$cell_id, meta$cell_id), all(meta$cell_type == "MPP"))
Y <- as.matrix(dat[, -1]); ids <- colnames(Y)
groups <- list(all_MPP = c(1, 2, 3), S1_S2_LK_branch = c(1, 2), S1_S3_CLP_branch = c(1, 3))
results <- list(); summaries <- list()

# Fit each MPP subset ----
for (group in names(groups)) {
  ix <- meta$State %in% groups[[group]]
  raw_t <- meta$Pseudotime[ix]
  t <- (raw_t - min(raw_t)) / diff(range(raw_t))
  batch <- droplevels(factor(meta$batch[ix]))
  A <- Y[ix, , drop = FALSE]
  grid <- seq(0, 1, length.out = 100)
  weights <- prop.table(table(batch))
  one <- function(j) {
    y <- A[, j]; warnings <- character()
    out <- list(feature_id = ids[j], group = group, n_cells = length(y), n_detected = sum(y > 0),
                fit_valid = FALSE, standardizable = FALSE, p_value = NA_real_,
                logLik_full = NA_real_, logLik_null = NA_real_, scale = NA_real_, reason = "",
                curve = rep(NA_real_, 100))
    tryCatch({
      if (any(tapply(y > 0, batch, sum) == 0)) stop("No detected values in at least one batch")
      d <- data.frame(y, t, batch)
      fit <- function(form) withCallingHandlers(
        survreg(form, data = d, dist = "gaussian", control = survreg.control(maxiter = 200)),
        warning = function(e) { warnings <<- c(warnings, conditionMessage(e)); invokeRestart("muffleWarning") })
      full <- fit(Surv(y, y > 0, type = "left") ~ batch + ns(t, df = 3))
      null <- fit(Surv(y, y > 0, type = "left") ~ batch)
      if (any(!is.finite(coef(full))) || any(!is.finite(coef(null))) ||
          !is.finite(full$scale) || full$scale <= 0 ||
          any(grepl("converge|infinite|singular", warnings, ignore.case = TRUE))) stop("Invalid or nonconvergent fit")
      lf <- as.numeric(logLik(full)); ln <- as.numeric(logLik(null))
      if (!is.finite(lf - ln) || lf < ln - 1e-5) stop("Invalid nested likelihood")
      mu <- rep(0, 100)
      for (b in levels(batch)) {
        g <- data.frame(t = grid, batch = factor(rep(b, 100), levels = levels(batch)))
        mu <- mu + as.numeric(weights[b]) * as.vector(predict(full, g, type = "lp"))
      }
      if (any(!is.finite(mu))) stop("Invalid prediction")
      out$fit_valid <- TRUE
      out$standardizable <- sd(mu) > 1e-8
      out$p_value <- pchisq(max(0, 2 * (lf - ln)), df = 3, lower.tail = FALSE)
      out$logLik_full <- lf; out$logLik_null <- ln; out$scale <- full$scale
      out$curve <- mu; out$reason <- paste(unique(warnings), collapse = "; ")
      out
    }, error = function(e) { out$reason <- conditionMessage(e); out })
  }
  results[[group]] <- lapply(seq_along(ids), one)
  write.csv(data.frame(cell_id = meta$cell_id[ix], State = meta$State[ix], batch = batch,
                       pseudotime = raw_t, normalized_pseudotime = t),
            file.path(output, paste0(group, "_cells.csv")), row.names = FALSE)
  summaries[[group]] <- data.frame(group = group, n_cells = sum(ix), n_S1 = sum(meta$State[ix] == 1),
                                   n_S2 = sum(meta$State[ix] == 2), n_S3 = sum(meta$State[ix] == 3),
                                   pseudotime_min = min(raw_t), pseudotime_max = max(raw_t))
  cat(group, "fitted:", length(ids), "features;", sum(ix), "cells\n")
}

# Adjust within each selected cell scope ----
stats <- do.call(rbind, lapply(results, function(rr) do.call(rbind, lapply(rr, function(x)
  as.data.frame(x[names(x) != "curve"], stringsAsFactors = FALSE)))))
stats$q_value <- NA_real_
ok <- stats$fit_valid & is.finite(stats$p_value)
for (group in names(groups)) {
  keep <- ok & stats$group == group
  stats$q_value[keep] <- p.adjust(stats$p_value[keep], method = "BH")
}
# Retain the across-scope correction for auditing; the viewer uses within-scope q.
stats$q_across_scopes <- NA_real_
stats$q_across_scopes[ok] <- p.adjust(stats$p_value[ok], method = "BH")
write.csv(stats, file.path(output, "all_fit_statistics.csv"), row.names = FALSE)
write.csv(do.call(rbind, summaries), file.path(output, "group_summary.csv"), row.names = FALSE)
for (group in names(results)) {
  mu <- do.call(rbind, lapply(results[[group]], `[[`, "curve"))
  colnames(mu) <- sprintf("t_%03d", 0:99)
  z <- t(apply(mu, 1, function(v) if (all(is.finite(v)) && sd(v) > 1e-8) (v - mean(v)) / sd(v) else rep(NA_real_, 100)))
  colnames(z) <- colnames(mu)
  write.csv(data.frame(feature_id = ids, mu, check.names = FALSE), file.path(output, paste0(group, "_predicted_log_curves.csv")), row.names = FALSE)
  write.csv(data.frame(feature_id = ids, z, check.names = FALSE), file.path(output, paste0(group, "_standardized_curves.csv")), row.names = FALSE)
}
print(with(stats, table(group, fit_valid)))
cat("Saved scope-specific fits and within-scope BH results.\n")
