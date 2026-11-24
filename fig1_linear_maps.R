# Fig 1: evolution under three linear genotype-phenotype maps x = B alpha.
#   A: single Gaussian peak, trajectories from four starting points.
#   B-D: two-peak Perlin landscape, outcomes from a grid of starting points
#        (10 diploid sexual replicates per point), basin border inferred by
#        SVM (red for B1, white for B2/B3), plus the Lande-equation separatrix
#        (dashed black) for comparison.

source("common.R")
suppressPackageStartupMessages({
  library(patchwork)
  library(ggnewscale)
})

B <- list(
  B1 = diag(2),
  B2 = matrix(c(5, 4, 4, 5) / 3, 2),
  B3 = matrix(c(5, -4, -4, 5) / 3, 2)
)
models <- lapply(B, gp_linear)

# Allelic values are unbounded for linear maps (B2 and B3 require negative
# alpha values in parts of phenotype space).
LOWER <- -Inf

# ---- A: single peak ----------------------------------------------------------
single <- landscape_gaussian(c(5.5, 5.5), omega = 3^2)
starts_A <- rbind(c(1, 5.5), c(5.5, 10), c(10, 10), c(10, 1))
trajA <- cached("fig1_single", run_jobs(nrow(starts_A) * 3L, function(k) {
  b <- (k - 1L) %/% nrow(starts_A) + 1L
  s <- (k - 1L) %% nrow(starts_A) + 1L
  args <- modifyList(PAPER, list(lower = LOWER, record_every = 10L))
  sim <- do.call(simulate_population, c(
    list(models[[b]], single, start = starts_A[s, ],
         stop = stop_at_peak_mean(matrix(c(5.5, 5.5), 1)), seed = 10L + k),
    args
  ))
  as_tibble(sim$trajectory) %>% mutate(map = names(B)[b], start = s)
}))

# ---- B-D: two-peak Perlin landscape ------------------------------------------
# Frequency-0.15 Perlin landscape (seed 96, offset (0, -2.5)) with two peaks on
# opposite sides of the starting grid x1 in [2, 8], x2 in [0, 6]:
# approximately (1.66, 2.93) and (9.34, 3.70).
FREQ <- 0.15
SEED_1 <- 96L
OFFSET <- c(0, -2.5)
perlin <- landscape_perlin(frequency = FREQ, seed = SEED_1, offset = OFFSET)
land <- landscape_grid(perlin, -10, 20, step = 0.02)
peaks <- find_peaks(perlin, -10, 20, res = 800)
in_view_peaks <- filter(peaks, between(x1, 0, 10), between(x2, 0, 7))
message(sprintf("Fig 1 two-peak landscape: frequency = %g, offset = (%g, %g), seed = %d, %d peaks in view",
                FREQ, OFFSET[1], OFFSET[2], SEED_1, nrow(in_view_peaks)))

GRID_N <- if (QUICK) 9L else 31L
REPS <- if (QUICK) 2L else 10L
grid <- expand.grid(x1 = seq(2, 8, length.out = GRID_N),
                    x2 = seq(0, 6, length.out = GRID_N))
jobs <- expand.grid(point = seq_len(nrow(grid)), map = seq_along(B), rep = seq_len(REPS))
res <- cached("fig1_two_peaks", run_jobs(nrow(jobs), function(k) {
  j <- jobs[k, ]
  run_to_peak(models[[j$map]], land, unlist(grid[j$point, ]), peaks,
              seed = 1000L + k, lower = LOWER) %>%
    mutate(map = names(B)[j$map], rep = j$rep)
}))
message(sprintf("Fig 1: %d runs, %d did not converge", nrow(res), sum(!res$converged)))

svm <- lapply(names(B), function(b) svm_borders(filter(res, map == b)))
names(svm) <- names(B)
borders <- bind_rows(lapply(names(B), function(b) mutate(svm[[b]]$borders, map = b)))
accuracy <- tibble(map = names(B), cv_accuracy = sapply(svm, `[[`, "accuracy"))
print(accuracy)
write_csv(accuracy, file.path(OUT_DIR, "fig1_svm_accuracy.csv"))

lande <- cached("fig1_lande", bind_rows(lapply(names(B), function(b) {
  mutate(lande_separatrices(models[[b]], perlin, lower = -2, upper = 12, step = 0.02),
         map = b)
}))) %>% filter(between(x1, 1.8, 8.2), between(x2, -0.2, 6.2))

# Exact diploid M-matrix entries for each linear map
M_mats <- lapply(models, function(m) {
  m_matrix(m, c(1, 1), mu = PAPER$mu, sigma = PAPER$sigma, ploidy = PAPER$ploidy)
})
M_header <- sapply(names(B), function(b) {
  M <- M_mats[[b]]
  sprintf("%s: M = [[%.6f, %.6f], [%.6f, %.6f]]",
          b, M[1, 1], M[1, 2], M[2, 1], M[2, 2])
})

# ---- Plotting ----------------------------------------------------------------
bgA <- landscape_df(single, 0, 11, res = 200)
map_cols <- c(B1 = "#F8766D", B2 = "#00BA38", B3 = "#619CFF")
map_labs <- c(B1 = expression(B[1] ~ "(unbiased)"),
              B2 = expression(B[2] ~ "(+ cov)"),
              B3 = expression(B[3] ~ "(- cov)"))

pA <- ggplot(bgA, aes(x, y)) +
  geom_raster(aes(fill = z)) +
  scale_fill_viridis_c(name = "Fitness", guide = "none") +
  geom_path(data = trajA,
            aes(mean_x1, mean_x2, colour = map, group = interaction(map, start)),
            linewidth = 0.85, inherit.aes = FALSE) +
  annotate("point", x = 5.5, y = 5.5, shape = 4, colour = "red", size = 2.5, stroke = 1) +
  scale_colour_manual(values = map_cols, labels = map_labs, name = "Linear map") +
  coord_equal(xlim = c(0, 11), ylim = c(0, 11), expand = FALSE) +
  labs(x = expression("Phenotype " * x[1]), y = expression("Phenotype " * x[2]))

agree <- res %>%
  count(map, start_x1, start_x2, peak) %>%
  group_by(map, start_x1, start_x2) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(freq = n / REPS)
ends <- res %>% mutate(px = peaks$x1[peak], py = peaks$x2[peak])
bg <- landscape_df(land, xlim = c(0, 10), ylim = c(-0.3, 7.0), res = 250)
b1_border <- filter(borders, map == "B1")

make_panel <- function(b, show_legend = FALSE) {
  M <- M_mats[[b]]
  p <- ggplot(bg, aes(x, y)) +
    geom_raster(aes(fill = z)) +
    scale_fill_viridis_c(name = "Fitness", limits = c(0, 1)) +
    geom_segment(data = filter(ends, map == b),
                 aes(x = start_x1, y = start_x2, xend = px, yend = py),
                 inherit.aes = FALSE, colour = "orange", linewidth = 0.12, alpha = 0.7) +
    geom_point(data = in_view_peaks, aes(x1, x2), colour = "red", shape = 4,
               size = 2.2, stroke = 1, inherit.aes = FALSE) +
    new_scale_colour() +
    geom_point(data = filter(agree, map == b),
               aes(start_x1, start_x2, colour = freq), size = 0.55, inherit.aes = FALSE) +
    scale_colour_gradient(low = "white", high = "black", limits = c(0.5, 1),
                          name = "Max freq.") +
    geom_tile(data = b1_border, aes(x1, x2), fill = "red",
              width = 0.06, height = 0.06, inherit.aes = FALSE)
  if (b != "B1") {
    p <- p + geom_tile(data = filter(borders, map == b), aes(x1, x2), fill = "white",
                       width = 0.06, height = 0.06, inherit.aes = FALSE)
  }
  p <- p +
    geom_path(data = filter(lande, map == b), aes(x1, x2, group = branch),
              colour = "black", linetype = 2, linewidth = 0.45, inherit.aes = FALSE) +
    coord_equal(xlim = c(0, 10), ylim = c(-0.3, 7.0), expand = FALSE) +
    labs(
      title = sprintf("%s: M = [%.6f, %.6f; %.6f, %.6f]",
                      b, M[1, 1], M[1, 2], M[2, 1], M[2, 2]),
      x = expression("Phenotype " * x[1]),
      y = expression("Phenotype " * x[2])
    ) +
    theme(plot.title = element_text(size = 9.5, face = "bold"))
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

pB <- make_panel("B1", show_legend = TRUE)
pC <- make_panel("B2", show_legend = FALSE)
pD <- make_panel("B3", show_legend = FALSE)

fig1 <- ((pA + labs(tag = "A")) | (pB + labs(tag = "B"))) /
        ((pC + labs(tag = "C")) | (pD + labs(tag = "D"))) +
  plot_layout(guides = "collect")

ggsave(file.path(FIG_DIR, "fig1.png"), fig1, width = 12, height = 8.5, dpi = 300, bg = "white")
message("Saved ", file.path(FIG_DIR, "fig1.png"))
