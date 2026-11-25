# Search for a robust Fig 6 example under diploid sexual populations:
# find a (landscape, start) where motifs 1 and 5 have nearly identical initial
# M-matrices, share the first segment of their evolutionary path, and then
# diverge reliably across replicates to different local peaks as M evolves.

source("common.R")

SEEDS <- seq_len(as.integer(Sys.getenv("N_SEEDS", if (QUICK) 10L else 80L)))
LIM <- 15
FAMILIES <- bind_rows(
  tibble(family = "perlin", frequency = c(0.1, 0.15, 0.2, 0.3)),
  tibble(family = "twin", frequency = c(0.05, 0.075, 0.1))
)
make_landscape <- function(family, frequency, seed) {
  if (family == "perlin") return(landscape_perlin(frequency = frequency, seed = seed, offset = 0.5))
  layer <- function(offset) {
    landscape_perlin(frequency = frequency, seed = seed, offset = offset,
                     transform = function(z) ((z + 1) / 2)^2)
  }
  landscape_smoothmax(layer(c(0, 1)), layer(c(-1, 0)), beta = 100)
}

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
Mf <- cached("fig6_search_mfields", list(`1` = m_field(1), `5` = m_field(5)))

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

shape <- function(M) {
  e <- eigen(M, symmetric = TRUE)
  c(angle = atan2(e$vectors[2, 1], e$vectors[1, 1]) %% pi,
    ecc = sqrt(1 - e$values[2] / e$values[1]),
    logdet = log(det(M)))
}
starts <- expand.grid(x1 = seq(0.5, 3, by = 0.5), x2 = seq(4, 12, by = 0.5))
sim_M <- starts %>%
  rowwise() %>%
  mutate(
    s1 = list(shape(M_at(1, c(x1, x2)))),
    s5 = list(shape(M_at(5, c(x1, x2)))),
    angle_diff = abs(((s1["angle"] - s5["angle"]) * 180 / pi + 90) %% 180 - 90),
    ecc1 = s1["ecc"], ecc5 = s5["ecc"],
    det_ratio = exp(s1["logdet"] - s5["logdet"])
  ) %>%
  ungroup() %>%
  select(-s1, -s5)
cand_starts <- sim_M %>% filter(angle_diff < 2, abs(ecc1 - ecc5) < 5e-3, between(det_ratio, 0.5, 2))
message(nrow(cand_starts), " candidate starting points with M1 ~ M5")

GSTEP <- 0.02
flow_path <- function(k, z, gx, gy, peaks, ds = 0.01, max_steps = 3000L) {
  path <- matrix(NA_real_, max_steps, 2)
  P <- as.matrix(peaks[, 1:2])
  dir <- function(z) {
    g <- c(bilin(gx, 0, GSTEP, z), bilin(gy, 0, GSTEP, z))
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
shared_length <- function(p1, p2, tol = 0.25) {
  n <- min(nrow(p1), nrow(p2))
  d <- sqrt(rowSums((p1[1:n, , drop = FALSE] - p2[1:n, , drop = FALSE])^2))
  split <- which(d > tol)
  (if (length(split)) split[1] else n) * 0.01
}

lands <- crossing(FAMILIES, seed = SEEDS)
flow <- cached("fig6_search_flow", run_jobs(nrow(lands), function(L) {
  l <- lands[L, ]
  perlin <- make_landscape(l$family, l$frequency, l$seed)
  xs <- seq(0, LIM, by = GSTEP)
  Z <- matrix(perlin(as.matrix(expand.grid(xs, xs))), length(xs))
  dr <- function(M) rbind(M[2, ] - M[1, ], (M[3:nrow(M), ] - M[1:(nrow(M) - 2), ]) / 2,
                          M[nrow(M), ] - M[nrow(M) - 1, ]) / GSTEP
  gx <- dr(Z); gy <- t(dr(t(Z)))
  peaks <- find_peaks(perlin, 0, LIM, res = 500)
  out <- list()
  for (r in seq_len(nrow(cand_starts))) {
    z0 <- c(cand_starts$x1[r], cand_starts$x2[r])
    f1 <- flow_path(1, z0, gx, gy, peaks)
    f5 <- flow_path(5, z0, gx, gy, peaks)
    if (is.na(f1$peak) || is.na(f5$peak) || f1$peak == f5$peak) next
    jit <- rbind(c(0.1, 0), c(-0.1, 0), c(0, 0.1), c(0, -0.1))
    robust <- sum(apply(jit, 1, function(d) {
      a <- flow_path(1, z0 + d, gx, gy, peaks)$peak
      b <- flow_path(5, z0 + d, gx, gy, peaks)$peak
      identical(a, f1$peak) && identical(b, f5$peak)
    })) >= 3L
    out[[length(out) + 1L]] <- tibble(
      family = l$family, frequency = l$frequency, seed = l$seed, x1 = z0[1], x2 = z0[2],
      peak1 = f1$peak, peak5 = f5$peak,
      peak1_x1 = peaks$x1[f1$peak], peak1_x2 = peaks$x2[f1$peak],
      peak5_x1 = peaks$x1[f5$peak], peak5_x2 = peaks$x2[f5$peak],
      peak_distance = sqrt((peaks$x1[f1$peak] - peaks$x1[f5$peak])^2 +
                           (peaks$x2[f1$peak] - peaks$x2[f5$peak])^2),
      shared = shared_length(f1$path, f5$path),
      len1 = nrow(f1$path) * 0.01, len5 = nrow(f5$path) * 0.01, robust = robust
    )
  }
  bind_rows(out)
})) %>% left_join(sim_M, by = c("x1", "x2"))

ranked <- flow %>%
  filter(robust, len1 < 15, len5 < 15) %>%
  mutate(score = shared / pmax(len1, len5)) %>%
  group_by(family, frequency, seed) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(score))
print(ranked %>% select(family, frequency, seed, x1, x2, angle_diff, det_ratio,
                        shared, len1, len5, peak_distance, score), n = 20)

TOP <- if (QUICK) 2L else 12L
REPS <- if (QUICK) 2L else 10L
top <- head(ranked, TOP)
jobs <- crossing(cand = seq_len(nrow(top)), motif = c(1L, 5L), rep = seq_len(REPS))
sims <- cached("fig6_search_sims", run_jobs(nrow(jobs), function(k) {
  j <- jobs[k, ]
  c0 <- top[j$cand, ]
  perlin <- make_landscape(c0$family, c0$frequency, c0$seed)
  land <- landscape_grid(perlin, 0, LIM)
  peaks <- find_peaks(perlin, 0, LIM, res = 500)
  run_to_peak(paper_motif(j$motif), land, c(c0$x1, c0$x2), peaks, seed = 7000L + k) %>%
    mutate(cand = j$cand, motif = j$motif, rep = j$rep,
           expected = if (j$motif == 1L) c0$peak1 else c0$peak5)
}))
summary_sims <- sims %>%
  group_by(cand, motif) %>%
  summarise(hit = mean(peak == expected), converged = mean(converged),
            generations = median(generations), .groups = "drop") %>%
  pivot_wider(names_from = motif, values_from = c(hit, converged, generations)) %>%
  left_join(mutate(top, cand = row_number()), by = "cand") %>%
  arrange(desc(pmin(hit_1, hit_5)), desc(score))
print(summary_sims %>% select(cand, family, frequency, seed, x1, x2,
                              hit_1, hit_5, generations_1, generations_5, shared, score))
write_csv(summary_sims, file.path(OUT_DIR, "fig6_search_candidates.csv"))
