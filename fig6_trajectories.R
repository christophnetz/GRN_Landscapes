# Fig 6: evolutionary trajectories of motifs 1 and 5 starting from nearly
# identical local M-matrices on a twin-peak adaptive landscape.
#
# Simulates diploid sexual populations (N = 10,000, split evenly between males
# and females) for 10,000 generations across REPS replicates per motif, tracking:
#   - mean phenotype (x1, x2),
#   - additive genetic covariance matrix G estimated via mid-parent-offspring
#     regression (random parent-offspring pairs independent of fitness),
#   - total phenotypic/genetic covariance matrix P,
#   - mutational covariance matrix M at the population mean genotype.

source("common.R")
suppressPackageStartupMessages(library(patchwork))

START <- c(1.0, 4.0)
REPS <- as.integer(Sys.getenv("FIG6_REPS", if (QUICK) 2L else 10L))
N_POP <- as.integer(Sys.getenv("FIG6_N", if (QUICK) 1000L else 10000L))
GENS <- as.integer(Sys.getenv("FIG6_GENS", if (QUICK) 2000L else 6000L))
REC_EVERY <- as.integer(Sys.getenv("FIG6_REC", 50L))
PO_PAIRS <- min(N_POP, 5000L)

# Two shifted, squared Perlin layers (frequency 0.05, seed 44) combined by a
# smooth maximum (beta = 100), selected via fig6_search.R for robust divergence
# of motifs 1 and 5 from nearly identical initial M-matrices.
layer <- function(offset) {
  landscape_perlin(frequency = 0.05, seed = 44L, offset = offset,
                   transform = function(z) ((z + 1) / 2)^2)
}
perlin <- landscape_smoothmax(layer(c(0, 1)), layer(c(-1, 0)), beta = 100)
land <- landscape_grid(perlin, 0, 15)
peaks <- find_peaks(perlin, 0, 15, res = 500)
in_view_peaks <- filter(peaks, between(x1, 0, 12), between(x2, 2, 12))

# Report initial M-matrix similarity at START
for (k in c(1L, 5L)) {
  Mk <- m_matrix(paper_motif(k), phenotype = START, mu = PAPER$mu,
                 sigma = PAPER$sigma, ploidy = PAPER$ploidy)
  ek <- eigen(Mk, symmetric = TRUE)
  ecc <- sqrt(1 - ek$values[2] / ek$values[1])
  ang <- atan2(ek$vectors[2, 1], ek$vectors[1, 1]) * 180 / pi
  message(sprintf("Initial M (motif %d) at (%.2f, %.2f): ecc=%.6f, angle=%.2f deg, det=%.3e, trace=%.3e",
                  k, START[1], START[2], ecc, ang, det(Mk), sum(diag(Mk))))
}

# Diploid sexual simulation matching grnland::simulate_population(), additionally
# recording parent-offspring regression G every record_every generations.
simulate_with_G_po <- function(model, landscape, start, N, generations,
                               mu, sigma, selection, record_every = 50L,
                               po_pairs = 5000L, seed = 1L) {
  set.seed(seed)
  L <- nrow(model$loci)
  traits <- model$traits
  g0 <- as.numeric(genotype_at(model, start))
  x_init <- model$x0; x_init[traits] <- start
  d0 <- develop(model, matrix(g0, 1), matrix(x_init, 1))
  if (d0$status != 0L) stop("starting genotype failed to develop")

  A <- array(rep(g0 / 2, each = N), c(N, L, 2L))
  X <- matrix(d0$X[1, ], N, model$n_genes, byrow = TRUE)
  W <- rep(landscape(X[1, traits, drop = FALSE]), N)
  bad <- rep(FALSE, N)
  fem <- seq_len(N) <= N %/% 2L

  n_rec <- generations %/% record_every + 2L
  rec <- matrix(NA_real_, n_rec, 15L, dimnames = list(NULL, c(
    "generation", "mean_x1", "mean_x2",
    "G_x1_x1", "G_x1_x2", "G_x2_x2",
    "P_x1_x1", "P_x1_x2", "P_x2_x2",
    "M_x1_x1", "M_x1_x2", "M_x2_x2",
    "mean_fitness", "n_invalid", "po_valid"
  )))
  idx_rec <- 0L

  record_state <- function(gen) {
    Xt <- X[!bad, traits, drop = FALSE]
    Pmat <- if (nrow(Xt) > 1L) stats::cov(Xt) else matrix(0, 2, 2)
    Gmat <- if (gen == 0L) matrix(0, 2, 2) else estimate_G_po(model, A, X, bad, n_pairs = po_pairs)
    gmean <- colMeans((A[, , 1] + A[, , 2])[!bad, , drop = FALSE])
    Mmat <- tryCatch(
      m_matrix(model, genotype = gmean, mu = mu, sigma = sigma, ploidy = 2L),
      error = function(e) matrix(NA_real_, 2, 2)
    )
    idx_rec <<- idx_rec + 1L
    rec[idx_rec, ] <<- c(
      gen, colMeans(Xt),
      Gmat[1, 1], Gmat[1, 2], Gmat[2, 2],
      Pmat[1, 1], Pmat[1, 2], Pmat[2, 2],
      Mmat[1, 1], Mmat[1, 2], Mmat[2, 2],
      mean(W, na.rm = TRUE), sum(bad), as.numeric(!anyNA(Gmat))
    )
  }

  record_state(0L)
  for (gen in seq_len(generations)) {
    w <- selection(W)
    mothers <- which(fem)[sample.int(sum(fem), N, TRUE, prob = w[fem])]
    fathers <- which(!fem)[sample.int(sum(!fem), N, TRUE, prob = w[!fem])]

    gamete <- function(par) {
      pick <- matrix(sample.int(2L, N * L, TRUE), N, L)
      matrix(A[cbind(rep(par, L), rep(seq_len(L), each = N), as.vector(pick))], N, L)
    }
    A <- array(c(gamete(mothers), gamete(fathers)), c(N, L, 2L))
    X <- X[mothers, , drop = FALSE]

    n_sites <- N * L * 2L
    nmut <- stats::rbinom(1L, n_sites, mu)
    if (nmut > 0L) {
      pos <- sample.int(n_sites, nmut)
      A[pos] <- pmax(A[pos] + stats::rnorm(nmut, 0, sigma), 0)
    }

    d <- develop(model, A[, , 1] + A[, , 2], X)
    ok <- d$status == 0L
    X <- d$X
    bad <- !ok
    W <- rep(NA_real_, N)
    if (any(ok)) W[ok] <- landscape(X[ok, traits, drop = FALSE])

    if (gen %% record_every == 0L) record_state(gen)
  }
  if (rec[idx_rec, "generation"] != generations) record_state(generations)
  as_tibble(rec[seq_len(idx_rec), , drop = FALSE])
}

jobs <- expand.grid(rep = seq_len(REPS), motif = c(1L, 5L))
res <- cached("fig6", run_jobs(nrow(jobs), function(k) {
  j <- jobs[k, ]
  tr <- simulate_with_G_po(
    paper_motif(j$motif), land, start = START, N = N_POP, generations = GENS,
    mu = PAPER$mu, sigma = PAPER$sigma, selection = PAPER$selection,
    record_every = REC_EVERY, po_pairs = PO_PAIRS, seed = 600L + k
  )
  last <- tr[nrow(tr), ]
  tr %>% mutate(motif = j$motif, rep = j$rep,
                final_peak = nearest_peak(c(last$mean_x1, last$mean_x2), peaks))
}, cores = min(N_CORES, nrow(jobs))))

outcomes <- res %>%
  distinct(motif, rep, final_peak) %>%
  mutate(peak_x1 = peaks$x1[final_peak], peak_x2 = peaks$x2[final_peak]) %>%
  count(motif, final_peak, peak_x1, peak_x2)
print(outcomes)
write_csv(outcomes, file.path(OUT_DIR, "fig6_outcomes.csv"))

# Pick a representative replicate per motif (preferring divergent pair if observed)
rep_peaks <- res %>% distinct(motif, rep, final_peak)
rep1 <- 1L
rep5 <- 1L
p1_main <- outcomes %>% filter(motif == 1) %>% slice_max(n, n = 1) %>% pull(final_peak)
p5_alt <- rep_peaks %>% filter(motif == 5, final_peak != p1_main[1])
if (nrow(p5_alt) > 0L) {
  rep1 <- rep_peaks %>% filter(motif == 1, final_peak == p1_main[1]) %>% pull(rep) %>% first()
  rep5 <- p5_alt$rep[1]
}
tr <- bind_rows(
  filter(res, motif == 1L, rep == rep1),
  filter(res, motif == 5L, rep == rep5)
)

# ---- A: landscape, all replicate trajectories, steepest ascent ---------------
ascent <- function(z, ds = 0.01, n = 3000L, h = 1e-4) {
  path <- matrix(NA_real_, n, 2)
  for (i in seq_len(n)) {
    path[i, ] <- z
    g <- (perlin(rbind(z + c(h, 0), z + c(0, h))) -
          perlin(rbind(z - c(h, 0), z - c(0, h)))) / (2 * h)
    if (sqrt(sum(g^2)) < 1e-6) break
    z <- z + ds * g / sqrt(sum(g^2))
  }
  as_tibble(path[complete.cases(path), , drop = FALSE], .name_repair = ~ c("x", "y"))
}
asc <- ascent(START)
bg <- landscape_df(land, xlim = c(0, 11.5), ylim = c(2.5, 11.5), res = 250)

pA <- ggplot(bg, aes(x, y)) +
  geom_raster(aes(fill = z)) +
  scale_fill_viridis_c(name = "Fitness") +
  geom_path(data = res, aes(mean_x1, mean_x2, group = rep),
            colour = "grey90", alpha = 0.5, linewidth = 0.35, inherit.aes = FALSE) +
  geom_path(data = asc, aes(x, y), colour = "red", linetype = 2, linewidth = 0.7,
            inherit.aes = FALSE) +
  geom_point(data = tr, aes(mean_x1, mean_x2, colour = generation),
             size = 0.7, inherit.aes = FALSE) +
  scale_colour_gradient(low = "#132B43", high = "#56B1F7", name = "Generation") +
  geom_point(data = in_view_peaks, aes(x1, x2), shape = 4, colour = "red",
             size = 2.2, stroke = 1, inherit.aes = FALSE) +
  scale_x_continuous(breaks = c(0, 3, 6, 9)) +
  scale_y_continuous(breaks = c(3, 6, 9)) +
  coord_equal(xlim = c(0, 11.5), ylim = c(2.5, 11.5), expand = FALSE) +
  facet_wrap(~motif, labeller = labeller(motif = motif_label)) +
  labs(x = expression("Phenotype " * x[1]), y = expression("Phenotype " * x[2])) +
  theme(panel.spacing.x = unit(0.9, "lines"))

# ---- B: trajectory with G (parent-offspring) and M (trace-matched) ellipses --
ell_step <- max(200L, REC_EVERY * 4L)
mats <- tr %>%
  filter(generation > 0L, generation %% ell_step == 0L) %>%
  mutate(
    G = pmap(list(G_x1_x1, G_x1_x2, G_x2_x2), ~ matrix(c(..1, ..2, ..2, ..3), 2)),
    M = pmap(list(M_x1_x1, M_x1_x2, M_x2_x2), ~ matrix(c(..1, ..2, ..2, ..3), 2)),
    M_scaled = map2(M, G, function(m, g) {
      if (anyNA(g) || sum(diag(g)) <= 0 || anyNA(m) || sum(diag(m)) <= 0) return(NULL)
      m * sum(diag(g)) / sum(diag(m))
    })
  )
build_ells <- function(col) {
  mats %>%
    mutate(e = map2(.data[[col]], map2(mean_x1, mean_x2, c), function(mat, ctr) {
      if (is.null(mat) || anyNA(mat)) return(NULL)
      ellipse_path(2^2 * mat, ctr)
    })) %>%
    select(motif, generation, e) %>%
    unnest(e)
}
pB <- ggplot(tr, aes(mean_x1, mean_x2)) +
  geom_path(linewidth = 0.5) +
  geom_path(data = build_ells("M_scaled"), aes(ex, ey, group = generation),
            colour = "#00BFC4", alpha = 0.65, linewidth = 0.45) +
  geom_path(data = build_ells("G"), aes(ex, ey, group = generation),
            colour = "#F8766D", alpha = 0.65, linewidth = 0.45) +
  geom_point(data = in_view_peaks, aes(x1, x2), shape = 4, colour = "red",
             size = 2.0, stroke = 0.9, inherit.aes = FALSE) +
  scale_x_continuous(breaks = c(0, 3, 6, 9)) +
  scale_y_continuous(breaks = c(3, 6, 9)) +
  coord_equal(xlim = c(0, 11.5), ylim = c(2.5, 11.5), expand = FALSE) +
  facet_wrap(~motif, labeller = labeller(motif = motif_label)) +
  labs(x = expression("Phenotype " * x[1]), y = expression("Phenotype " * x[2])) +
  theme(panel.spacing.x = unit(0.9, "lines"))

# ---- C: eccentricity, angle, determinant of G and M --------------------------
summ <- tr %>%
  filter(generation > 0L) %>%
  pivot_longer(c(G_x1_x1, G_x1_x2, G_x2_x2, M_x1_x1, M_x1_x2, M_x2_x2),
               names_to = c("matrix", ".value"), names_pattern = "([GM])_(.*)") %>%
  mutate(
    matrix = if_else(matrix == "G", "G-matrix", "M-matrix"),
    disc = sqrt((x1_x1 - x2_x2)^2 + 4 * x1_x2^2),
    lmax = (x1_x1 + x2_x2 + disc) / 2,
    lmin = (x1_x1 + x2_x2 - disc) / 2,
    eccentricity = if_else(lmax > 0, sqrt(pmax(0, 1 - pmax(lmin, 0) / lmax)), NA_real_),
    angle = 0.5 * atan2(2 * x1_x2, x1_x1 - x2_x2) * 180 / pi,
    determinant = pmax(x1_x1 * x2_x2 - x1_x2^2, 1e-12),
    total_var = pmax(x1_x1 + x2_x2, 1e-12)
  )

make_c <- function(var, ylab, log_y = FALSE, show_strip = FALSE, show_x = FALSE) {
  p <- ggplot(summ, aes(generation, .data[[var]], colour = matrix)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = c("G-matrix" = "#F8766D", "M-matrix" = "#00BFC4"),
                        name = NULL) +
    facet_wrap(~motif, labeller = labeller(motif = motif_label)) +
    labs(x = if (show_x) "Generation" else NULL, y = ylab) +
    theme(panel.spacing.x = unit(1.0, "lines"))
  if (log_y) p <- p + scale_y_log10()
  if (!show_strip) p <- p + theme(strip.text = element_blank())
  if (!show_x) p <- p + theme(axis.text.x = element_blank())
  p
}

pC <- ((make_c("eccentricity", "Eccentricity", show_strip = TRUE) + labs(tag = "C")) /
       make_c("angle", "Angle (deg)") /
       make_c("determinant", "Determinant", log_y = TRUE, show_x = TRUE)) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

fig6 <- ((pA + labs(tag = "A")) / (pB + labs(tag = "B")) | pC) +
  plot_layout(widths = c(1.25, 1))

ggsave(file.path(FIG_DIR, "fig6.png"), fig6, width = 12.5, height = 8.0, dpi = 300, bg = "white")
message("Saved ", file.path(FIG_DIR, "fig6.png"))
