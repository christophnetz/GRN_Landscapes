# Fig 4 and Fig 5: effect of adaptive landscape complexity on evolutionary
# outcomes and alignment with developmental bias (M_max).
#
# For each Perlin frequency in FREQS, simulate N_LAND random landscapes; each
# of the 6 motifs is run once from every point on a grid over [2, 10]^2.
#   Fig 4A: mean outcome mismatch with motif 6 (no developmental bias) vs frequency.
#   Fig 4B: proportion of identical outcomes between all pairs of motifs.
#   Fig 5:  angle between evolutionary divergence vector (start -> peak reached)
#           and dominant eigenvector M_max of M at the reached peak.
#
# Environment variables:
#   SIGMA=0.2   mutational sd (default 0.2 as in Figs 4-5; set SIGMA=0.05 for
#               the sensitivity check reported in the text)

source("common.R")
suppressPackageStartupMessages(library(patchwork))

SIGMA <- as.numeric(Sys.getenv("SIGMA", "0.2"))
FREQS <- if (QUICK) c(0.3, 0.6, 1.8) else c(0.3, 0.6, 1.2, 1.8)
N_LAND <- as.integer(Sys.getenv("FIG4_N_LAND", if (QUICK) 2L else 20L))
GRID_N <- as.integer(Sys.getenv("FIG4_GRID_N", if (QUICK) 9L else 25L))
MAX_GENS <- as.integer(Sys.getenv("FIG4_MAX_GENS", 3000L))

set.seed(20261117L)
lands <- crossing(frequency = FREQS, landscape = seq_len(N_LAND)) %>%
  mutate(seed = sample.int(100000L, n()))
grid <- expand.grid(x1 = seq(2, 10, length.out = GRID_N),
                    x2 = seq(2, 10, length.out = GRID_N))
tag <- paste0("sigma", SIGMA)

make_land <- function(l) {
  perlin <- landscape_perlin(frequency = l$frequency, seed = l$seed, offset = 0.5)
  list(perlin = perlin,
       land = landscape_grid(perlin, 0, 30, step = 0.02),
       peaks = find_peaks(perlin, 0, 30, res = 1200))
}

res <- cached(paste0("fig4_", tag), {
  bind_rows(lapply(seq_len(nrow(lands)), function(L) {
    l <- lands[L, ]
    message(sprintf("Fig 4/5 (sigma=%g): frequency %g, landscape %d/%d (seed %d)",
                    SIGMA, l$frequency, l$landscape, N_LAND, l$seed))
    ld <- make_land(l)
    n_peaks_box <- sum(between(ld$peaks$x1, 1, 11) & between(ld$peaks$x2, 1, 11))
    jobs <- expand.grid(point = seq_len(nrow(grid)), motif = 1:6)
    run_jobs(nrow(jobs), function(k) {
      j <- jobs[k, ]
      run_to_peak(paper_motif(j$motif), ld$land, unlist(grid[j$point, ]), ld$peaks,
                  seed = 1000000L * L + k, sigma = SIGMA, generations = MAX_GENS) %>%
        mutate(motif = j$motif)
    }) %>%
      mutate(frequency = l$frequency, landscape = l$landscape,
             n_peaks_box = n_peaks_box,
             peak_x1 = ld$peaks$x1[peak], peak_x2 = ld$peaks$x2[peak])
  }))
})
message(sprintf("Fig 4/5 (%s): %d runs, %d did not converge (%.2f%%)",
                tag, nrow(res), sum(!res$converged), 100 * mean(!res$converged)))

# ---- Fig 4: outcome mismatch -------------------------------------------------
# Retain only starting points where all 6 motifs converged.
wide <- res %>%
  group_by(frequency, landscape, start_x1, start_x2) %>%
  filter(all(converged)) %>%
  ungroup() %>%
  select(frequency, landscape, n_peaks_box, start_x1, start_x2, motif, peak, peak_x1, peak_x2) %>%
  pivot_wider(names_from = motif, values_from = c(peak, peak_x1, peak_x2), names_sep = "_m")

pairs <- combn(1:6, 2, simplify = FALSE)
pairwise <- bind_rows(lapply(pairs, function(p) {
  c1 <- paste0("peak_m", p[1]); c2 <- paste0("peak_m", p[2])
  x1a <- paste0("peak_x1_m", p[1]); y1a <- paste0("peak_x2_m", p[1])
  x1b <- paste0("peak_x1_m", p[2]); y1b <- paste0("peak_x2_m", p[2])
  wide %>%
    group_by(frequency, landscape, n_peaks_box) %>%
    summarise(
      same = mean(.data[[c1]] == .data[[c2]]),
      mean_dist = mean(sqrt((.data[[x1a]] - .data[[x1b]])^2 +
                            (.data[[y1a]] - .data[[y1b]])^2)),
      n_points = n(),
      .groups = "drop"
    ) %>%
    mutate(motif_1 = p[1], motif_2 = p[2])
}))
write_csv(pairwise, file.path(OUT_DIR, paste0("fig4_pairwise_", tag, ".csv")))

# Summary table for manuscript reporting
summary_mismatch <- pairwise %>%
  filter(motif_2 == 6) %>%
  group_by(frequency, motif_1) %>%
  summarise(mean_mismatch = mean(1 - same), sd_mismatch = sd(1 - same),
            mean_dist = mean(mean_dist), .groups = "drop")
print(summary_mismatch, n = 30)
write_csv(summary_mismatch, file.path(OUT_DIR, paste0("fig4_summary_", tag, ".csv")))

motif_cols <- c("1" = "#F8766D", "2" = "#A3A500", "3" = "#00BF7D",
                "4" = "#00B0F6", "5" = "#E76BF3", "6" = "#FF61CC")

pA <- pairwise %>%
  filter(motif_2 == 6) %>%
  mutate(motif = factor(motif_1)) %>%
  ggplot(aes(factor(frequency), 1 - same, colour = motif)) +
  geom_point(position = position_dodge(width = 0.55), alpha = 0.75, size = 1.6) +
  stat_summary(aes(group = motif), fun = mean, geom = "point", shape = 18,
               size = 3.2, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = motif_cols[1:5], labels = motif_label, name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Landscape complexity", y = "Outcome mismatch with motif 6")

pooled <- pairwise %>%
  group_by(motif_1, motif_2) %>%
  summarise(same = mean(same), .groups = "drop")

pB <- if (requireNamespace("ggraph", quietly = TRUE)) {
  g <- igraph::graph_from_data_frame(
    pooled %>% transmute(from = motif_label(motif_1), to = motif_label(motif_2), same),
    directed = FALSE
  )
  set.seed(1)
  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(aes(width = same, alpha = same)) +
    ggraph::geom_node_point(size = 5) +
    ggraph::geom_node_label(aes(label = name), repel = TRUE, size = 3.8,
                            fill = "white", alpha = 0.85, label.size = 0) +
    ggraph::scale_edge_width(range = c(0.4, 3.2), name = "Same outcomes") +
    ggraph::scale_edge_alpha(range = c(0.25, 1), name = "Same outcomes") +
    theme_void()
} else {
  ggplot(pooled, aes(factor(motif_1), factor(motif_2), fill = same)) +
    geom_tile() +
    labs(x = "Motif", y = "Motif", fill = "Same outcomes")
}

fig4 <- (pA + labs(tag = "A")) + (pB + labs(tag = "B")) + plot_layout(widths = c(1.2, 1))
ggsave(file.path(FIG_DIR, paste0("fig4_", tag, ".png")), fig4,
       width = 11.5, height = 5.2, dpi = 300, bg = "white")
if (SIGMA == 0.2) {
  ggsave(file.path(FIG_DIR, "fig4.png"), fig4, width = 11.5, height = 5.2, dpi = 300, bg = "white")
}

# ---- Fig 5: alignment of evolutionary divergence with M_max ------------------
# Divergence vector: start -> peak reached; compared with dominant eigenvector
# of M at the evolutionary end point (peak reached).
peak_M <- res %>%
  filter(converged) %>%
  distinct(motif, peak_x1, peak_x2)
peak_M$angle_vec <- mapply(function(k, x1, x2) {
  M <- tryCatch(
    m_matrix(paper_motif(k), phenotype = c(x1, x2), mu = PAPER$mu,
             sigma = SIGMA, ploidy = PAPER$ploidy),
    error = function(e) NULL
  )
  if (is.null(M)) return(NA_real_)
  v <- eigen(M, symmetric = TRUE)$vectors[, 1]
  atan2(v[2], v[1])
}, peak_M$motif, peak_M$peak_x1, peak_M$peak_x2)

# Note: for motif 6, M is exactly isotropic, so eigen() returns an arbitrary
# coordinate axis; the uniform angle distribution reflects the isotropic
# distribution of divergence directions rather than lack of alignment per se.
align <- res %>%
  filter(converged) %>%
  left_join(peak_M, by = c("motif", "peak_x1", "peak_x2")) %>%
  mutate(dx = peak_x1 - start_x1, dy = peak_x2 - start_x2) %>%
  filter(sqrt(dx^2 + dy^2) > 0.1, !is.na(angle_vec)) %>%
  mutate(angle = acos(pmin(1, abs(cos(atan2(dy, dx) - angle_vec)))) * 180 / pi)

summary_align <- align %>%
  group_by(frequency, motif) %>%
  summarise(median_angle = median(angle), iqr_angle = IQR(angle),
            n = n(), .groups = "drop")
print(summary_align, n = 30)
write_csv(summary_align, file.path(OUT_DIR, paste0("fig5_summary_", tag, ".csv")))

p5 <- ggplot(align, aes(factor(frequency), angle, fill = factor(motif))) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.25, linewidth = 0.35,
               position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 45, linetype = 2, colour = "red", linewidth = 0.6) +
  scale_fill_manual(values = motif_cols, name = "Motif") +
  scale_y_continuous(breaks = seq(0, 90, by = 15),
                     labels = function(x) paste0(x, "\u00b0")) +
  labs(x = "Landscape complexity",
       y = expression("Angle between divergence and " * M[max]))

ggsave(file.path(FIG_DIR, paste0("fig5_", tag, ".png")), p5,
       width = 9.5, height = 5.5, dpi = 300, bg = "white")
if (SIGMA == 0.2) {
  ggsave(file.path(FIG_DIR, "fig5.png"), p5, width = 9.5, height = 5.5, dpi = 300, bg = "white")
}

# ---- Fig 5 alternative: relative mutational variance along divergence -------
# Alternative to the M_max angle (from the original boxplotsfig5.r):
#   e_d = (d' M d / d' d) / (tr(M) / 2),
# the mutational variance in the direction of divergence d (start -> peak
# reached) relative to the mean mutational variance over all directions.
# e_d > 1: evolution proceeded along a direction of above-average mutational
# variance; e_d < 1: below average. Isotropic M (motif 6) gives exactly 1, and
# e_d is in [lambda_min, lambda_max] / mean(lambda), so it weights alignment by
# the eccentricity of M. Computed with M at the end point (peak reached) and
# at the starting point.
m_entries <- function(k, x1, x2) {
  M <- tryCatch(
    m_matrix(paper_motif(k), phenotype = c(x1, x2), mu = PAPER$mu,
             sigma = SIGMA, ploidy = PAPER$ploidy),
    error = function(e) matrix(NA_real_, 2, 2)
  )
  c(M[1, 1], M[1, 2], M[2, 2])
}
add_M <- function(d, x1, x2, prefix) {
  E <- t(mapply(m_entries, d$motif, d[[x1]], d[[x2]]))
  d[[paste0(prefix, "11")]] <- E[, 1]
  d[[paste0(prefix, "12")]] <- E[, 2]
  d[[paste0(prefix, "22")]] <- E[, 3]
  d
}
rel_var <- function(dx, dy, m11, m12, m22) {
  (dx^2 * m11 + 2 * dx * dy * m12 + dy^2 * m22) / (dx^2 + dy^2) /
    ((m11 + m22) / 2)
}
M_end <- add_M(distinct(align, motif, peak_x1, peak_x2), "peak_x1", "peak_x2", "Me")
M_start <- add_M(distinct(align, motif, start_x1, start_x2), "start_x1", "start_x2", "Ms")
e_d <- align %>%
  left_join(M_end, by = c("motif", "peak_x1", "peak_x2")) %>%
  left_join(M_start, by = c("motif", "start_x1", "start_x2")) %>%
  transmute(frequency, landscape, motif, start_x1, start_x2, peak_x1, peak_x2,
            `M at end point` = rel_var(dx, dy, Me11, Me12, Me22),
            `M at starting point` = rel_var(dx, dy, Ms11, Ms12, Ms22)) %>%
  pivot_longer(c(`M at end point`, `M at starting point`),
               names_to = "M_at", values_to = "e_d") %>%
  filter(is.finite(e_d))

summary_e_d <- e_d %>%
  group_by(M_at, frequency, motif) %>%
  summarise(median_e_d = median(e_d), iqr_e_d = IQR(e_d),
            frac_above_1 = mean(e_d > 1 + 1e-8), n = n(), .groups = "drop")
print(summary_e_d, n = 60)
write_csv(summary_e_d, file.path(OUT_DIR, paste0("fig5_e_d_summary_", tag, ".csv")))

p5_ed <- ggplot(e_d, aes(factor(frequency), e_d, fill = factor(motif))) +
  geom_hline(yintercept = 1, linetype = 2, colour = "red", linewidth = 0.6) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.25, linewidth = 0.35,
               position = position_dodge(width = 0.8)) +
  facet_wrap(~M_at, ncol = 1) +
  scale_fill_manual(values = motif_cols, name = "Motif") +
  labs(x = "Landscape complexity",
       y = "Relative mutational variance\nin the direction of divergence")

ggsave(file.path(FIG_DIR, paste0("fig5_e_d_", tag, ".png")), p5_ed,
       width = 9.5, height = 8, dpi = 300, bg = "white")
message("Saved Fig 4 and Fig 5 for ", tag)
