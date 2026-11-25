# Shared settings and helper functions for all manuscript figure scripts.
# Run from the repository root with BLAS threads capped, e.g.:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript fig2_m_matrices.R
#
# Environment variables:
#   QUICK=1               reduced grids and replicates for fast pipeline checks
#   REDO=1                ignore cached simulation output and recompute
#   N_CORES=k             number of parallel workers
#   SELECTION=relative_min|proportional
#                         selection scheme (default: relative_min)

suppressPackageStartupMessages({
  library(grnland)
  library(tidyverse)
  library(parallel)
})

QUICK <- Sys.getenv("QUICK") == "1"
REDO <- Sys.getenv("REDO") == "1"
N_CORES <- as.integer(Sys.getenv("N_CORES", min(max(1L, detectCores() - 1L), 128L)))
SELECTION_NAME <- Sys.getenv("SELECTION", "relative_min")

SELECTION_FUN <- switch(
  SELECTION_NAME,
  relative_min = selection_relative_min(0.001),
  proportional = selection_proportional(),
  stop("Unknown SELECTION: ", SELECTION_NAME)
)

sub_tag <- c(
  if (QUICK) "quick",
  if (SELECTION_NAME != "relative_min") SELECTION_NAME
)
OUT_DIR <- if (length(sub_tag)) file.path("output", paste(sub_tag, collapse = "_")) else "output"
FIG_DIR <- if (length(sub_tag)) file.path("figures", paste(sub_tag, collapse = "_")) else "figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Simulation settings matching the manuscript: diploid sexual populations
# split evenly between males and females, two additive loci (alpha_1, alpha_2),
# per-allele mutation rate mu = 0.025, mutational effect sd sigma = 0.05.
PAPER <- list(
  N = 1000L,
  generations = 10000L,
  mu = 0.025,
  sigma = 0.05,
  ploidy = 2L,
  reproduction = "sexual",
  sexes = TRUE,
  selection = SELECTION_FUN
)

# Run fun(k) for k in seq_len(n) in parallel. Each job sets its own seed
# inside fun, so results do not depend on the number of cores.
run_jobs <- function(n, fun, cores = N_CORES) {
  cores <- max(1L, min(as.integer(cores), as.integer(n)))
  res <- mclapply(seq_len(n), function(k) {
    tryCatch(fun(k), error = function(e) {
      message("job ", k, " failed: ", conditionMessage(e))
      NULL
    })
  }, mc.cores = cores, mc.preschedule = FALSE)
  bind_rows(res)
}

# Load cached output or compute and save it.
cached <- function(name, expr) {
  path <- file.path(OUT_DIR, paste0(name, ".rds"))
  if (file.exists(path) && !REDO) return(readRDS(path))
  t0 <- Sys.time()
  val <- expr
  message(sprintf("%s computed in %.2f min", name,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(val, path)
  val
}

# Stop when the population mean phenotype remains within Euclidean distance
# `radius` of the same peak for `holds` consecutive checkpoints (`every`
# generations apart). Works cleanly across both small (sigma = 0.05) and large
# (sigma = 0.2) mutational scales and under weaker proportional selection.
stop_at_peak_mean <- function(peaks, radius = 0.1, holds = 2L, every = 10L) {
  P <- as.matrix(peaks)[, 1:2, drop = FALSE]
  last_k <- NA_integer_
  streak <- 0L
  r2 <- radius^2
  fun <- function(X) {
    if (nrow(X) == 0L) return(NULL)
    m <- colMeans(X)
    d2 <- colSums((t(P) - m)^2)
    k <- which.min(d2)
    if (d2[k] <= r2) {
      if (identical(k, last_k)) streak <<- streak + 1L else { last_k <<- k; streak <<- 1L }
      if (streak >= holds) return(k)
    } else {
      last_k <<- NA_integer_
      streak <<- 0L
    }
    NULL
  }
  structure(list(fun = fun, every = every, peaks = P), class = "grn_stop")
}

# One evolutionary run from a monomorphic start; returns the peak reached.
# Runs that do not meet the stopping rule within the generation limit are
# assigned the peak nearest to their final mean and flagged with converged = FALSE.
run_to_peak <- function(model, land, start, peaks, seed, ...) {
  args <- modifyList(PAPER, list(...))
  sim <- do.call(simulate_population, c(
    list(model, land, start = start,
         stop = stop_at_peak_mean(peaks), seed = seed,
         record_every = args$generations),
    args
  ))
  last <- sim$trajectory[nrow(sim$trajectory), ]
  end <- c(last$mean_x1, last$mean_x2)
  tibble(
    start_x1 = start[1], start_x2 = start[2],
    end_x1 = end[1], end_x2 = end[2],
    peak = if (is.na(sim$outcome)) nearest_peak(end, peaks) else sim$outcome,
    converged = !is.na(sim$outcome),
    generations = sim$generations
  )
}

# Basin borders from a radial SVM fitted to (start, outcome), plus the
# 10-fold cross-validated classification accuracy (%).
svm_borders <- function(d, res = 200, cost = 10, gamma = 0.5) {
  d <- d %>% transmute(x1 = start_x1, x2 = start_x2, peak = factor(peak))
  if (n_distinct(d$peak) < 2) {
    return(list(borders = tibble(x1 = numeric(), x2 = numeric()), accuracy = 100))
  }
  set.seed(1)
  fit <- e1071::svm(peak ~ x1 + x2, data = d, kernel = "radial",
                    cost = cost, gamma = gamma, cross = 10)
  xs <- seq(min(d$x1), max(d$x1), length.out = res)
  ys <- seq(min(d$x2), max(d$x2), length.out = res)
  g <- expand.grid(x1 = xs, x2 = ys)
  P <- matrix(as.integer(predict(fit, g)), res, res)
  B <- matrix(FALSE, res, res)
  B[-res, ] <- B[-res, ] | P[-1, ] != P[-res, ]
  B[, -res] <- B[, -res] | P[, -1] != P[, -res]
  idx <- which(B, arr.ind = TRUE)
  list(borders = tibble(x1 = xs[idx[, 1]], x2 = ys[idx[, 2]]),
       accuracy = fit$tot.accuracy)
}

# Landscape values on a grid, for plotting.
landscape_df <- function(land, lower = 0, upper = 11.5, res = 300,
                         xlim = c(lower, upper), ylim = c(lower, upper)) {
  g <- expand.grid(
    x = seq(xlim[1], xlim[2], length.out = res),
    y = seq(ylim[1], ylim[2], length.out = res)
  )
  g$z <- land(as.matrix(g))
  g
}

# Covariance ellipse (level-quantile contour of a bivariate normal).
ellipse_path <- function(M, centre, level = 0.95, n = 60) {
  if (any(!is.finite(M))) return(NULL)
  e <- eigen(M, symmetric = TRUE)
  if (any(e$values <= 0)) return(NULL)
  t <- seq(0, 2 * pi, length.out = n)
  r <- sqrt(qchisq(level, 2))
  pts <- r * cbind(cos(t), sin(t)) %*% diag(sqrt(e$values)) %*% t(e$vectors)
  tibble(ex = pts[, 1] + centre[1], ey = pts[, 2] + centre[2])
}

# Estimate additive genetic covariance matrix (G) via mid-parent-offspring
# covariance: independently sample random mother-father pairs irrespective of
# fitness, generate one diploid offspring per pair with free recombination and
# no mutation, develop offspring, and compute G = Cov(z_O, z_mid) + Cov(z_O, z_mid)^T
# (since Cov(z_O, z_mid) = 0.5 * G in expectation).
estimate_G_po <- function(model, A, X, bad, n_pairs = NULL) {
  N <- nrow(X)
  L <- dim(A)[2]
  traits <- model$traits
  nt <- length(traits)
  if (is.null(n_pairs)) n_pairs <- N
  fem <- which(!bad & seq_len(N) <= N %/% 2L)
  mal <- which(!bad & seq_len(N) > N %/% 2L)
  if (length(fem) < 2L || length(mal) < 2L) return(matrix(NA_real_, nt, nt))

  mothers <- sample(fem, n_pairs, replace = TRUE)
  fathers <- sample(mal, n_pairs, replace = TRUE)

  gamete <- function(par) {
    pick <- matrix(sample.int(2L, n_pairs * L, replace = TRUE), n_pairs, L)
    matrix(A[cbind(rep(par, L), rep(seq_len(L), each = n_pairs), as.vector(pick))],
           n_pairs, L)
  }
  G_off <- gamete(mothers) + gamete(fathers)
  d_off <- develop(model, G_off, X[mothers, , drop = FALSE])
  ok <- d_off$status == 0L
  if (sum(ok) < 3L) return(matrix(NA_real_, nt, nt))

  Z_mid <- 0.5 * (X[mothers[ok], traits, drop = FALSE] +
                  X[fathers[ok], traits, drop = FALSE])
  Z_off <- d_off$X[ok, traits, drop = FALSE]
  C_op <- stats::cov(Z_off, Z_mid)
  C_op + t(C_op)
}

motif_label <- function(x) paste("Motif", x)

theme_set(theme_minimal())
