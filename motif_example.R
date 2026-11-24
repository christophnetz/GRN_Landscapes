# Motif divergence example: all six motifs started from the same phenotype on
# the same rugged landscape evolve to different peaks because their M-matrices
# differ in shape and orientation.
#
# Step 1 (deterministic): screen Perlin landscapes and starting points under the
#   Lande flow dz/dt = M_k(z) grad W(z), keeping robust starts where the six
#   motifs reach >= 3 distinct peaks.
# Step 2 (stochastic): simulate top candidates with diploid sexual replicate
#   populations.
# Step 3: plot the best example (default: freq = 0.4, seed = 6, start = (4, 8)).

source("common.R")
suppressPackageStartupMessages(library(patchwork))

LIM <- 20
MOTIFS <- 1:6
FREQS <- c(0.3, 0.4)
SEEDS <- seq_len(as.integer(Sys.getenv("N_SEEDS", if (QUICK) 6L else 40L)))
TOP <- if (QUICK) 2L else 8L
REPS <- if (QUICK) 3L else 10L

# ---- M fields of all motifs, precomputed on a grid --------------------------
MSTEP <- 0.05
mx <- seq(MSTEP, LIM, by = MSTEP)
m_field <- function(k) {
  m <- paper_motif(k)
  g <- expand.grid(x1 = mx, x2 = mx)
  vals <- t(mapply(function(a, b) {
    M <- tryCatch(
      m_matrix(m, phenotype = c(a, b), mu = PAPER$mu, sigma = PAPER$sigma,
               ploidy = PAPER$ploidy),
      error = function(e) matrix(NA, 2, 2)
    )
    c(M[1, 1], M[1, 2], M[2, 2])
  }, g$x1, g$x2))
  lapply(1:3, function(i) matrix(vals[, i], length(mx)))
}
Mf <- cached("motif_example_mfields",
             setNames(mclapply(MOTIFS, m_field, mc.cores = min(N_CORES, 6L)), MOTIFS))

bilin <- function(Z, lo, step, z) {
  n1 <- nrow(Z); n2 <- ncol(Z)
  u <- min(max((z[1] - lo) / step, 0), n1 - 1)
  v <- min(max((z[2] - lo) / step, 0), n2 - 1)
  i <- min(floor(u), n1 - 2); j <- min(floor(v), n2 - 2)
  wu <- u - i; wv <- v - j
  (1 - wu) * (1 - wv) * Z[i + 1, j + 1] + wu * (1 - wv) * Z[i + 2, j + 1] +
    (1 - wu) * wv * Z[i + 1, j + 2] + wu * wv * Z[i + 2, j + 2]
}
M_at <- function(k, z) {
  f <- Mf[[as.character(k)]]
  a <- bilin(f[[1]], MSTEP, MSTEP, z)
  b <- bilin(f[[2]], MSTEP, MSTEP, z)
  d <- bilin(f[[3]], MSTEP, MSTEP, z)
  matrix(c(a, b, b, d), 2)
}

# ---- Lande flow --------------------------------------------------------------
GSTEP <- 0.02
gradient_grid <- function(land) {
  xs <- seq(0, LIM, by = GSTEP)
  Z <- matrix(land(as.matrix(expand.grid(xs, xs))), length(xs))
  dr <- function(M) rbind(M[2, ] - M[1, ], (M[3:nrow(M), ] - M[1:(nrow(M) - 2), ]) / 2,
                          M[nrow(M), ] - M[nrow(M) - 1, ]) / GSTEP
  list(gx = dr(Z), gy = t(dr(t(Z))))
}
flow_path <- function(k, z, grads, peaks, ds = 0.01, max_steps = 2000L) {
  path <- matrix(NA_real_, max_steps, 2)
  P <- as.matrix(peaks[, 1:2])
  dir <- function(z) {
    g <- c(bilin(grads$gx, 0, GSTEP, z), bilin(grads$gy, 0, GSTEP, z))
    v <- as.numeric(M_at(k, z) %*% g)
    v / sqrt(sum(v^2))
  }
  for (i in seq_len(max_steps)) {
    path[i, ] <- z
    d2 <- colSums((t(P) - z)^2)
    if (min(d2) < 0.03^2) return(list(path = path[1:i, , drop = FALSE], peak = which.min(d2)))
    k1 <- dir(z)
    k2 <- dir(z + 0.5 * ds * k1)
    if (any(!is.finite(k2))) break
    z <- z + ds * k2
    if (any(z < MSTEP) || any(z > LIM)) break
  }
  list(path = path[complete.cases(path), , drop = FALSE], peak = NA_integer_)
}
flow_peaks <- function(z, grads, peaks) vapply(MOTIFS, function(k) flow_path(k, z, grads, peaks)$peak, 1L)

# ---- Step 1: deterministic screen --------------------------------------------
starts <- expand.grid(x1 = 2:12, x2 = 2:12)
lands <- crossing(frequency = FREQS, seed = SEEDS)
JIT <- rbind(c(0.15, 0), c(-0.15, 0), c(0, 0.15), c(0, -0.15))
screen <- cached("motif_example_flow", run_jobs(nrow(lands), function(L) {
  l <- lands[L, ]
  land <- landscape_perlin(frequency = l$frequency, seed = l$seed)
  grads <- gradient_grid(land)
  peaks <- find_peaks(land, 0, LIM, res = 500)
  out <- list()
  for (r in seq_len(nrow(starts))) {
    z0 <- c(starts$x1[r], starts$x2[r])
    pk <- flow_peaks(z0, grads, peaks)
    if (anyNA(pk) || n_distinct(pk) < 3L) next
    jp <- t(apply(JIT, 1, function(d) flow_peaks(z0 + d, grads, peaks)))
    agree <- colMeans(jp == matrix(pk, nrow(JIT), length(MOTIFS), byrow = TRUE), na.rm = TRUE)
    out[[length(out) + 1L]] <- tibble(
      frequency = l$frequency, seed = l$seed, x1 = z0[1], x2 = z0[2],
      n_peaks = n_distinct(pk), min_agree = min(agree),
      peaks = list(pk), peak_x1 = list(peaks$x1[pk]), peak_x2 = list(peaks$x2[pk])
    )
  }
  bind_rows(out)
}))

ranked <- screen %>%
  filter(min_agree >= 0.75,
         map_lgl(peak_x1, ~ all(.x < 15)), map_lgl(peak_x2, ~ all(.x < 15))) %>%
  group_by(frequency, seed) %>%
  arrange(desc(n_peaks), desc(min_agree), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(desc(n_peaks), desc(min_agree))

# ---- Step 2: stochastic diploid sexual check ---------------------------------
top <- head(ranked, TOP)
jobs <- crossing(cand = seq_len(nrow(top)), motif = MOTIFS, rep = seq_len(REPS))
sims <- cached("motif_example_sims", run_jobs(nrow(jobs), function(k) {
  j <- jobs[k, ]
  c0 <- top[j$cand, ]
  perlin <- landscape_perlin(frequency = c0$frequency, seed = c0$seed)
  peaks <- find_peaks(perlin, 0, LIM, res = 500)
  args <- modifyList(PAPER, list(record_every = 20L))
  sim <- do.call(simulate_population, c(
    list(paper_motif(j$motif), landscape_grid(perlin, 0, LIM),
         start = c(c0$x1, c0$x2), stop = stop_at_peak_mean(peaks), seed = 8000L + k),
    args
  ))
  last <- sim$trajectory[nrow(sim$trajectory), ]
  as_tibble(sim$trajectory) %>%
    select(generation, mean_x1, mean_x2) %>%
    mutate(
      cand = j$cand, motif = j$motif, rep = j$rep,
      peak = if (is.na(sim$outcome)) nearest_peak(c(last$mean_x1, last$mean_x2), peaks) else sim$outcome,
      converged = !is.na(sim$outcome),
      predicted = c0$peaks[[1]][j$motif]
    )
}))

outcomes <- sims %>% distinct(cand, motif, rep, peak, converged, predicted)
summary_sims <- outcomes %>%
  group_by(cand, motif) %>%
  summarise(hit = mean(peak == predicted), converged = mean(converged), .groups = "drop") %>%
  group_by(cand) %>%
  summarise(min_hit = min(hit), mean_hit = mean(hit), converged = mean(converged), .groups = "drop") %>%
  left_join(mutate(top, cand = row_number()), by = "cand") %>%
  arrange(desc(n_peaks * (min_hit >= 0.6)), desc(mean_hit))
print(summary_sims %>% select(cand, frequency, seed, x1, x2, n_peaks, min_hit, mean_hit, converged))
write_csv(select(summary_sims, -peaks, -peak_x1, -peak_x2),
          file.path(OUT_DIR, "motif_example_candidates.csv"))

# ---- Step 3: figure ----------------------------------------------------------
motif_names <- c(
  `1` = "1: x1 \u2192 x2", `2` = "2: x1 \u22a3 x2",
  `3` = "3: x1 \u2192 x2 \u22a3 x1", `4` = "4: x1 \u22a3 x2 \u22a3 x1",
  `5` = "5: x1 \u2192 x2 \u2192 x1", `6` = "6: none"
)
cols <- setNames(c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", "#D55E00"), MOTIFS)

plot_example <- function(cand_row, file) {
  perlin <- landscape_perlin(frequency = cand_row$frequency, seed = cand_row$seed)
  peaks <- find_peaks(perlin, 0, LIM, res = 500)
  start <- c(cand_row$x1, cand_row$x2)
  grads <- gradient_grid(perlin)
  lande <- bind_rows(lapply(MOTIFS, function(k) {
    p <- flow_path(k, start, grads, peaks)$path
    tibble(x1 = p[, 1], x2 = p[, 2], motif = k)
  }))
  tr <- filter(sims, cand == cand_row$cand)
  out <- filter(outcomes, cand == cand_row$cand)
  reached <- sort(unique(out$peak))
  peak_lab <- tibble(peak = reached, label = LETTERS[seq_along(reached)],
                     x1 = peaks$x1[reached], x2 = peaks$x2[reached])

  xs <- c(start[1], peak_lab$x1, tr$mean_x1, lande$x1)
  ys <- c(start[2], peak_lab$x2, tr$mean_x2, lande$x2)
  half <- max(diff(range(xs)), diff(range(ys))) / 2 + 1.2
  cx <- mean(range(xs)); cy <- mean(range(ys))
  xlim <- c(cx - half, cx + half); ylim <- c(cy - half, cy + half)
  bg <- expand.grid(x = seq(xlim[1], xlim[2], length.out = 300),
                    y = seq(ylim[1], ylim[2], length.out = 300))
  bg$z <- perlin(as.matrix(bg))

  pA <- ggplot(bg, aes(x, y)) +
    geom_raster(aes(fill = z)) +
    scale_fill_gradient(low = "grey10", high = "grey90", name = "Fitness") +
    geom_path(data = tr, aes(mean_x1, mean_x2, colour = factor(motif),
                             group = interaction(motif, rep)),
              linewidth = 0.3, alpha = 0.6, inherit.aes = FALSE) +
    geom_path(data = lande, aes(x1, x2, colour = factor(motif)),
              linewidth = 1, linetype = "42", inherit.aes = FALSE) +
    geom_point(data = filter(peaks, between(x1, xlim[1], xlim[2]),
                             between(x2, ylim[1], ylim[2])),
               aes(x1, x2), shape = 4, size = 2, stroke = 1, colour = "white",
               inherit.aes = FALSE) +
    geom_label(data = peak_lab, aes(x1, x2, label = label), nudge_x = 0.45,
               nudge_y = 0.45, size = 3.5, label.padding = unit(0.12, "lines"),
               inherit.aes = FALSE) +
    annotate("point", start[1], start[2], shape = 21, size = 3.5,
             fill = "white", colour = "black", stroke = 1) +
    scale_colour_manual(values = cols, labels = motif_names, name = "Motif") +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = expression(x[1]), y = expression(x[2]),
         subtitle = sprintf("Perlin frequency %.1f, seed %d; start (%g, %g)\nthin: %d diploid sexual replicates per motif; dashed: Lande flow",
                            cand_row$frequency, cand_row$seed, start[1], start[2], REPS))

  ell <- bind_rows(lapply(MOTIFS, function(k) {
    M <- m_matrix(paper_motif(k), phenotype = start, mu = PAPER$mu,
                  sigma = PAPER$sigma, ploidy = PAPER$ploidy)
    ellipse_path(M / sum(diag(M)), c(0, 0), level = 0.5) %>%
      mutate(motif = k, label = motif_names[as.character(k)])
  }))
  pB <- ggplot(ell, aes(ex, ey)) +
    geom_polygon(aes(fill = factor(motif)), colour = "black", linewidth = 0.3) +
    scale_fill_manual(values = cols, guide = "none") +
    facet_wrap(~label, nrow = 2) +
    coord_equal() +
    labs(x = NULL, y = NULL, subtitle = "M at the start, scaled to equal total variance") +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          strip.text = element_text(size = 7.0),
          panel.background = element_rect(fill = "grey95", colour = NA))

  pC <- out %>%
    left_join(peak_lab, by = "peak") %>%
    ggplot(aes(factor(motif), fill = label)) +
    geom_bar() +
    scale_fill_grey(start = 0.25, end = 0.85) +
    labs(x = "Motif", y = "Replicates", fill = "Peak", subtitle = "Peak reached")

  ggsave(file, (pA + labs(tag = "A")) |
           ((pB + labs(tag = "B")) / (pC + labs(tag = "C")) + plot_layout(heights = c(3, 2))),
         width = 13, height = 7, dpi = 300, bg = "white")
}

perfect <- filter(summary_sims, min_hit == 1)
if (nrow(perfect) == 0L) perfect <- summary_sims[1, ]
default <- filter(perfect, frequency == 0.4, seed == 6)
pick <- as.integer(Sys.getenv("EXAMPLE", if (nrow(default)) default$cand[1] else perfect$cand[1]))
best <- filter(summary_sims, cand == pick)
message(sprintf("Example: frequency %.1f, seed %d, start (%g, %g), min_hit = %.2f",
                best$frequency, best$seed, best$x1, best$x2, best$min_hit))
plot_example(best, file.path(FIG_DIR, "motif_example.png"))
