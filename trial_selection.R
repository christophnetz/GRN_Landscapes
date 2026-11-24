# Trial comparing selection schemes under diploid sexual populations:
#   - relative_min:  w - min(w) + 0.001 (strong range-dependent selection)
#   - proportional:  weights proportional to w (standard Wright-Fisher)
#
# Evaluates both schemes on:
#   1) Complex Perlin landscape (freq 0.5, seed 1245), motifs 1, 2, 6
#   2) Complexity sweep (freqs 0.3, 0.6, 1.2, 1.8; 5 landscapes; sigma = 0.2)
#   3) Motif divergence example (freq 0.4, seed 6, start (4, 8))

source("common.R")
suppressPackageStartupMessages(library(patchwork))

SCHEMES <- list(
  relative_min = selection_relative_min(0.001),
  proportional = selection_proportional()
)

# ---- 1. Complex landscape basins (Motifs 1, 2, 6) ----------------------------
perlin_c <- landscape_perlin(frequency = 0.5, seed = 1245L)
land_c <- landscape_grid(perlin_c, 0, 15)
peaks_c <- find_peaks(perlin_c, 0, 15, res = 500)
GRID_N <- if (QUICK) 9L else 21L
REPS_B <- if (QUICK) 2L else 5L
grid_b <- expand.grid(x1 = seq(2, 10, length.out = GRID_N),
                      x2 = seq(2, 10, length.out = GRID_N))

trial_basins <- cached("trial_selection_basins", {
  jobs <- expand.grid(point = seq_len(nrow(grid_b)), motif = c(1L, 2L, 6L),
                      rep = seq_len(REPS_B), scheme = names(SCHEMES),
                      stringsAsFactors = FALSE)
  run_jobs(nrow(jobs), function(k) {
    j <- jobs[k, ]
    run_to_peak(paper_motif(j$motif), land_c, unlist(grid_b[j$point, ]), peaks_c,
                seed = 50000L + k, selection = SCHEMES[[j$scheme]]) %>%
      mutate(motif = j$motif, rep = j$rep, scheme = j$scheme)
  })
})

basin_summary <- trial_basins %>%
  group_by(scheme, motif) %>%
  summarise(
    pct_converged = 100 * mean(converged),
    median_gen = median(generations),
    .groups = "drop"
  )
print(basin_summary)

# Compare outcomes between the two schemes at identical starting points
basin_agree <- trial_basins %>%
  group_by(scheme, motif, start_x1, start_x2) %>%
  count(peak) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(scheme, motif, start_x1, start_x2, peak) %>%
  pivot_wider(names_from = scheme, values_from = peak) %>%
  group_by(motif) %>%
  summarise(scheme_agreement = mean(relative_min == proportional, na.rm = TRUE),
            .groups = "drop")
print(basin_agree)

# ---- 2. Complexity sweep (outcome mismatch & M_max alignment) ----------------
FREQS <- c(0.3, 0.6, 1.2, 1.8)
N_LAND_T <- if (QUICK) 2L else 5L
GRID_T <- if (QUICK) 9L else 17L
set.seed(20261117L)
lands_t <- crossing(frequency = FREQS, landscape = seq_len(N_LAND_T)) %>%
  mutate(seed = sample.int(100000L, n()))
grid_t <- expand.grid(x1 = seq(2, 10, length.out = GRID_T),
                      x2 = seq(2, 10, length.out = GRID_T))

trial_comp <- cached("trial_selection_complexity", {
  bind_rows(lapply(seq_len(nrow(lands_t)), function(L) {
    l <- lands_t[L, ]
    perlin <- landscape_perlin(frequency = l$frequency, seed = l$seed, offset = 0.5)
    land <- landscape_grid(perlin, 0, 30, step = 0.02)
    peaks <- find_peaks(perlin, 0, 30, res = 1200)
    jobs <- expand.grid(point = seq_len(nrow(grid_t)), motif = 1:6,
                        scheme = names(SCHEMES), stringsAsFactors = FALSE)
    run_jobs(nrow(jobs), function(k) {
      j <- jobs[k, ]
      run_to_peak(paper_motif(j$motif), land, unlist(grid_t[j$point, ]), peaks,
                  seed = 900000L * L + k, sigma = 0.2, generations = 3000L,
                  selection = SCHEMES[[j$scheme]]) %>%
        mutate(motif = j$motif, scheme = j$scheme)
    }) %>%
      mutate(frequency = l$frequency, landscape = l$landscape,
             peak_x1 = peaks$x1[peak], peak_x2 = peaks$x2[peak])
  }))
})

comp_conv <- trial_comp %>%
  group_by(scheme, frequency) %>%
  summarise(pct_converged = 100 * mean(converged),
            median_gen = median(generations), .groups = "drop")
print(comp_conv)

wide_t <- trial_comp %>%
  group_by(scheme, frequency, landscape, start_x1, start_x2) %>%
  filter(all(converged)) %>%
  ungroup() %>%
  select(scheme, frequency, landscape, start_x1, start_x2, motif, peak) %>%
  pivot_wider(names_from = motif, values_from = peak, names_prefix = "m")

mm_t <- bind_rows(lapply(1:5, function(k) {
  wide_t %>%
    group_by(scheme, frequency, landscape) %>%
    summarise(mismatch = mean(.data[[paste0("m", k)]] != m6), .groups = "drop") %>%
    mutate(motif = k)
}))

mm_summary <- mm_t %>%
  group_by(scheme, frequency, motif) %>%
  summarise(mean_mismatch = mean(mismatch), .groups = "drop")
print(mm_summary, n = 40)
write_csv(mm_summary, file.path(OUT_DIR, "trial_selection_mismatch.csv"))

# Alignment angles
peak_M_t <- trial_comp %>%
  filter(converged) %>%
  distinct(motif, peak_x1, peak_x2)
peak_M_t$angle_vec <- mapply(function(k, x1, x2) {
  M <- tryCatch(m_matrix(paper_motif(k), phenotype = c(x1, x2), mu = PAPER$mu,
                         sigma = 0.2, ploidy = PAPER$ploidy),
                error = function(e) NULL)
  if (is.null(M)) return(NA_real_)
  v <- eigen(M, symmetric = TRUE)$vectors[, 1]
  atan2(v[2], v[1])
}, peak_M_t$motif, peak_M_t$peak_x1, peak_M_t$peak_x2)

align_t <- trial_comp %>%
  filter(converged) %>%
  left_join(peak_M_t, by = c("motif", "peak_x1", "peak_x2")) %>%
  mutate(dx = peak_x1 - start_x1, dy = peak_x2 - start_x2) %>%
  filter(sqrt(dx^2 + dy^2) > 0.1, !is.na(angle_vec)) %>%
  mutate(angle = acos(pmin(1, abs(cos(atan2(dy, dx) - angle_vec)))) * 180 / pi)

align_summary <- align_t %>%
  group_by(scheme, frequency, motif) %>%
  summarise(median_angle = median(angle), .groups = "drop")
print(align_summary, n = 50)
write_csv(align_summary, file.path(OUT_DIR, "trial_selection_alignment.csv"))

# ---- Plot comparison ---------------------------------------------------------
p1 <- ggplot(mm_summary, aes(factor(frequency), mean_mismatch,
                             colour = factor(motif), group = motif)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~scheme) +
  labs(x = "Landscape complexity", y = "Mean outcome mismatch with motif 6",
       colour = "Motif", title = "Outcome mismatch: relative_min vs proportional")

p2 <- ggplot(align_summary, aes(factor(frequency), median_angle,
                                colour = factor(motif), group = motif)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 45, linetype = 2, colour = "grey40") +
  facet_wrap(~scheme) +
  labs(x = "Landscape complexity", y = expression("Median angle to " * M[max] * " (deg)"),
       colour = "Motif", title = "Alignment with M_max: relative_min vs proportional")

ggsave(file.path(FIG_DIR, "trial_selection_comparison.png"), p1 / p2,
       width = 10, height = 7.5, dpi = 300, bg = "white")
message("Saved ", file.path(FIG_DIR, "trial_selection_comparison.png"))
