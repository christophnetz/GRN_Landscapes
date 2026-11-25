# Fig 2: Gene regulatory motifs, M-matrices, and sensitivity vectors over
# phenotype space. Computes the diploid mutational variance-covariance matrix
# M = 2 * mu * sigma^2 * S S' (matching the simulation settings) and magnifies
# the 95% confidence ellipse axes by an explicit factor MAG.

source("common.R")
suppressPackageStartupMessages(library(patchwork))

MAG <- 26   # magnification of 95% ellipse axes (standard deviations)
pts <- expand.grid(x1 = seq(2, 10, length.out = 5),
                   x2 = seq(2, 10, length.out = 5))

rows <- list()
for (k in 1:6) {
  m <- paper_motif(k)
  for (r in seq_len(nrow(pts))) {
    x <- unlist(pts[r, ])
    g <- genotype_at(m, x)
    S <- sensitivity(m, g)$S
    M <- m_matrix(m, g, mu = PAPER$mu, sigma = PAPER$sigma, ploidy = PAPER$ploidy)
    rows[[length(rows) + 1L]] <- tibble(
      motif = k, x1 = x[1], x2 = x[2],
      s1_x = S[1, 1], s1_y = S[2, 1], s2_x = S[1, 2], s2_y = S[2, 2],
      ell = list(ellipse_path(MAG^2 * M, x))
    )
  }
}
d <- bind_rows(rows)
ell <- d %>% select(motif, x1, x2, ell) %>% unnest(ell)

# Pure-R schematic diagram of each 2-gene regulatory motif
draw_motif_diagram <- function(k) {
  circle <- function(cx, cy, r = 0.55, n = 80) {
    t <- seq(0, 2 * pi, length.out = n)
    tibble(x = cx + r * cos(t), y = cy + r * sin(t))
  }
  c1 <- circle(-1.3, 0)
  c2 <- circle(1.3, 0)

  p <- ggplot() +
    geom_polygon(data = c1, aes(x, y), fill = "white", colour = "black", linewidth = 0.5) +
    geom_polygon(data = c2, aes(x, y), fill = "white", colour = "black", linewidth = 0.5) +
    annotate("text", x = -1.3, y = 0, label = "x[1]", parse = TRUE, size = 3.5) +
    annotate("text", x = 1.3, y = 0, label = "x[2]", parse = TRUE, size = 3.5) +
    coord_equal(xlim = c(-2.3, 2.3), ylim = c(-0.85, 0.85), expand = FALSE) +
    labs(title = paste("Motif", k)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10.5,
                                    margin = margin(b = 2)))

  # Helper: draw activation arrow or repression T-bar from x_start to x_end at height y
  add_edge <- function(p, x0, x1, y, type = c("act", "rep")) {
    type <- match.arg(type)
    if (type == "act") {
      p + annotate("segment", x = x0, xend = x1, y = y, yend = y,
                   linewidth = 0.55,
                   arrow = arrow(length = unit(0.14, "cm"), type = "open"))
    } else {
      p + annotate("segment", x = x0, xend = x1, y = y, yend = y, linewidth = 0.55) +
        annotate("segment", x = x1, xend = x1, y = y - 0.14, yend = y + 0.14,
                 linewidth = 0.65)
    }
  }

  if (k == 1) {
    p <- add_edge(p, -0.62, 0.62, 0, "act")
  } else if (k == 2) {
    p <- add_edge(p, -0.62, 0.62, 0, "rep")
  } else if (k == 3) {
    p <- add_edge(p, -0.60, 0.60, 0.18, "act")
    p <- add_edge(p, 0.60, -0.60, -0.18, "rep")
  } else if (k == 4) {
    p <- add_edge(p, -0.60, 0.60, 0.18, "rep")
    p <- add_edge(p, 0.60, -0.60, -0.18, "rep")
  } else if (k == 5) {
    p <- add_edge(p, -0.60, 0.60, 0.18, "act")
    p <- add_edge(p, 0.60, -0.60, -0.18, "act")
  }
  p
}

xlims <- range(c(ell$ex, d$x1 + d$s1_x, d$x1 + d$s2_x))
ylims <- range(c(ell$ey, d$x2 + d$s1_y, d$x2 + d$s2_y))

draw_motif_field <- function(k, show_x = FALSE, show_y = FALSE) {
  dk <- filter(d, motif == k)
  ek <- filter(ell, motif == k)
  ggplot(ek, aes(ex, ey, group = interaction(x1, x2))) +
    geom_path(linewidth = 0.45) +
    geom_segment(data = dk, aes(x = x1, y = x2, xend = x1 + s1_x, yend = x2 + s1_y),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
                 colour = "#2166AC", linewidth = 0.45) +
    geom_segment(data = dk, aes(x = x1, y = x2, xend = x1 + s2_x, yend = x2 + s2_y),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
                 colour = "#D6604D", linewidth = 0.45) +
    coord_equal(xlim = xlims, ylim = ylims) +
    scale_x_continuous(breaks = c(5, 10)) +
    scale_y_continuous(breaks = c(5, 10)) +
    labs(x = if (show_x) expression(x[1]) else NULL,
         y = if (show_y) expression(x[2]) else NULL) +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = if (show_y) element_text() else element_blank(),
          axis.text.x = if (show_x) element_text() else element_blank())
}

panels <- lapply(1:6, function(k) {
  diag_p <- draw_motif_diagram(k)
  field_p <- draw_motif_field(k, show_x = (k >= 4), show_y = (k %in% c(1, 4)))
  diag_p / field_p + plot_layout(heights = c(0.24, 1))
})

fig2 <- (panels[[1]] | panels[[2]] | panels[[3]]) /
        (panels[[4]] | panels[[5]] | panels[[6]])

ggsave(file.path(FIG_DIR, "fig2.png"), fig2, width = 9.5, height = 7.5, dpi = 300, bg = "white")
message(sprintf("Fig 2: 95%% ellipses of M (mu = %g, sigma = %g, ploidy = %d), axes magnified %gx",
                PAPER$mu, PAPER$sigma, PAPER$ploidy, MAG))
