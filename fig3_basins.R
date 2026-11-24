# Fig 3, Fig S1, Fig S2: evolutionary outcomes from a grid of starting points
# for all six motifs on a simple (Perlin frequency 0.1, seed 1234) and a
# complex (frequency 0.5, seed 1245) adaptive landscape, basin borders inferred
# by SVM (white; motif 6 reference in red), and separatrices predicted by the
# Lande equation with G proportional to M (black in Figs S1 and S2).

source("common.R")
suppressPackageStartupMessages({
  library(patchwork)
  library(ggnewscale)
})

LANDSCAPES <- tibble(
  name = c("simple", "complex"),
  label = c("Landscape complexity 0.1", "Landscape complexity 0.5"),
  frequency = c(0.1, 0.5),
  seed = c(1234L, 1245L)
)
GRID_N <- as.integer(Sys.getenv("FIG3_GRID_N", if (QUICK) 9L else 41L))
REPS <- as.integer(Sys.getenv("FIG3_REPS", if (QUICK) 2L else 10L))
grid <- expand.grid(x1 = seq(2, 10, length.out = GRID_N),
                    x2 = seq(2, 10, length.out = GRID_N))

fig3_rows <- list()

for (L in seq_len(nrow(LANDSCAPES))) {
  ls <- LANDSCAPES[L, ]
  perlin <- landscape_perlin(frequency = ls$frequency, seed = ls$seed)
  land <- landscape_grid(perlin, -5, 25, step = 0.02)
  peaks <- find_peaks(perlin, -5, 25, res = 800)

  jobs <- expand.grid(point = seq_len(nrow(grid)), motif = 1:6, rep = seq_len(REPS))
  res <- cached(paste0("fig3_", ls$name), run_jobs(nrow(jobs), function(k) {
    j <- jobs[k, ]
    run_to_peak(paper_motif(j$motif), land, unlist(grid[j$point, ]), peaks,
                seed = 100000L * L + k) %>%
      mutate(motif = j$motif, rep = j$rep)
  }))
  message(sprintf("Fig 3 (%s): %d runs, %d did not converge (%.2f%%)",
                  ls$name, nrow(res), sum(!res$converged), 100 * mean(!res$converged)))

  svm <- lapply(1:6, function(k) svm_borders(filter(res, motif == k)))
  borders <- bind_rows(lapply(1:6, function(k) mutate(svm[[k]]$borders, motif = k)))
  accuracy <- tibble(landscape = ls$name, motif = 1:6,
                     cv_accuracy = sapply(svm, `[[`, "accuracy"))
  print(accuracy)
  write_csv(accuracy, file.path(OUT_DIR, paste0("fig3_svm_accuracy_", ls$name, ".csv")))

  lande <- cached(paste0("fig3_lande_", ls$name), bind_rows(mclapply(1:6, function(k) {
    mutate(lande_separatrices(paper_motif(k), perlin, lower = 0, upper = 12, step = 0.01),
           motif = k)
  }, mc.cores = min(6L, N_CORES)))) %>%
    filter(between(x1, 2, 10), between(x2, 2, 10))

  agree <- res %>%
    count(motif, start_x1, start_x2, peak) %>%
    group_by(motif, start_x1, start_x2) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(freq = n / REPS)
  ends <- res %>% mutate(px = peaks$x1[peak], py = peaks$x2[peak])
  bg <- landscape_df(land, 0, 12, res = 300)

  outcome_plot <- function(motifs, ref_border = TRUE, lande_lines = FALSE,
                           row_strip = NULL, show_x = TRUE) {
    bg_p <- bg
    ends_p <- filter(ends, motif %in% motifs)
    pts_p <- ends_p %>% distinct(motif, px, py)
    agree_p <- filter(agree, motif %in% motifs)
    bord_p <- filter(borders, motif %in% motifs)
    lande_p <- filter(lande, motif %in% motifs)
    if (!is.null(row_strip)) {
      bg_p$landscape <- row_strip
      ends_p$landscape <- row_strip
      pts_p$landscape <- row_strip
      agree_p$landscape <- row_strip
      bord_p$landscape <- row_strip
      lande_p$landscape <- row_strip
    }
    p <- ggplot(bg_p, aes(x, y)) +
      geom_raster(aes(fill = z)) +
      scale_fill_viridis_c(name = "Fitness", limits = c(0, 1)) +
      geom_segment(data = ends_p,
                   aes(x = start_x1, y = start_x2, xend = px, yend = py),
                   inherit.aes = FALSE, colour = "orange", linewidth = 0.1, alpha = 0.7) +
      geom_point(data = pts_p,
                 aes(px, py), inherit.aes = FALSE, shape = 4, size = 1.6,
                 stroke = 0.9, colour = "red")
    if (ref_border) {
      ref <- filter(borders, motif == 6) %>% select(-motif) %>% crossing(motif = motifs)
      if (!is.null(row_strip)) ref$landscape <- row_strip
      p <- p + geom_tile(data = ref, aes(x1, x2), fill = "red",
                         width = 0.06, height = 0.06, inherit.aes = FALSE)
    }
    p <- p +
      new_scale_colour() +
      geom_point(data = agree_p,
                 aes(start_x1, start_x2, colour = freq), size = 0.35, inherit.aes = FALSE) +
      scale_colour_gradient(low = "white", high = "black", limits = c(0.5, 1),
                            name = "Max freq.") +
      geom_tile(data = bord_p, aes(x1, x2), fill = "white",
                width = 0.06, height = 0.06, inherit.aes = FALSE)
    if (lande_lines) {
      p <- p + geom_path(data = lande_p,
                         aes(x1, x2, group = branch), colour = "black",
                         linewidth = 0.35, inherit.aes = FALSE)
    }
    p <- p +
      coord_equal(xlim = c(0, 12), ylim = c(0, 12), expand = FALSE) +
      scale_x_continuous(breaks = c(0, 3, 6, 9, 12)) +
      scale_y_continuous(breaks = c(0, 3, 6, 9, 12)) +
      labs(x = if (show_x) expression("Phenotype " * x[1]) else NULL,
           y = expression("Phenotype " * x[2])) +
      theme(panel.spacing.x = unit(0.9, "lines"))
    if (!is.null(row_strip)) {
      p <- p + facet_grid(landscape ~ motif,
                          labeller = labeller(motif = motif_label))
    } else {
      p <- p + facet_wrap(~motif, labeller = labeller(motif = motif_label))
    }
    p
  }

  fig3_rows[[L]] <- outcome_plot(c(1, 2), ref_border = TRUE, lande_lines = FALSE,
                                 row_strip = ls$label, show_x = (L == nrow(LANDSCAPES)))

  # Individual row images (for vignette compatibility) and full S1/S2 grids
  ggsave(file.path(FIG_DIR, paste0("fig3_", ls$name, ".png")),
         outcome_plot(c(1, 2), ref_border = TRUE, lande_lines = FALSE),
         width = 8.5, height = 4.2, dpi = 300, bg = "white")
  ggsave(file.path(FIG_DIR, paste0("figS", L, ".png")),
         outcome_plot(1:6, ref_border = FALSE, lande_lines = TRUE),
         width = 11, height = 7.5, dpi = 300, bg = "white")
}

fig3_combined <- (fig3_rows[[1]] / fig3_rows[[2]]) +
  plot_layout(guides = "collect")
ggsave(file.path(FIG_DIR, "fig3.png"), fig3_combined, width = 9.5, height = 8.5, dpi = 300, bg = "white")
message("Saved ", file.path(FIG_DIR, "fig3.png"), ", figS1.png, figS2.png")
