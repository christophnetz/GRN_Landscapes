# Fig S3: Perlin adaptive landscapes of increasing frequency (identical seed)
# and their spatial autocorrelation profiles.

source("common.R")
suppressPackageStartupMessages(library(patchwork))

FREQS <- c(0.1, 0.3, 0.6, 1.2, 1.8)
SEED <- 1234L

maps <- bind_rows(lapply(FREQS, function(f) {
  landscape_df(landscape_perlin(frequency = f, seed = SEED, offset = 0.5),
               0, 15, res = 250) %>%
    mutate(frequency = f)
}))

pA <- ggplot(maps, aes(x, y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Fitness", limits = c(0, 1)) +
  coord_equal(expand = FALSE) +
  facet_wrap(~frequency, nrow = 1,
             labeller = labeller(frequency = function(f) paste("Frequency:", f))) +
  labs(x = expression(x[1]), y = expression(x[2])) +
  theme(panel.spacing.x = unit(0.8, "lines"))

# Spatial autocorrelation: correlation of fitness between pairs of points at
# distance d in random directions across a large region.
set.seed(1)
n <- 25000L
dists <- seq(0, 15, by = 0.15)
acf <- cached("figS3_acf", bind_rows(lapply(FREQS, function(f) {
  land <- landscape_perlin(frequency = f, seed = SEED, offset = 0.5)
  X <- matrix(runif(2L * n, 0, 50), ncol = 2L)
  w0 <- land(X)
  tibble(frequency = f, distance = dists, correlation = sapply(dists, function(d) {
    a <- runif(n, 0, 2 * pi)
    cor(w0, land(X + d * cbind(cos(a), sin(a))))
  }))
})))

pB <- ggplot(acf, aes(distance, correlation, colour = factor(frequency))) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  labs(x = "Distance", y = "Autocorrelation", colour = "Frequency")

figS3 <- (pA + labs(tag = "A")) / (pB + labs(tag = "B")) +
  plot_layout(heights = c(1, 0.9))

ggsave(file.path(FIG_DIR, "figS3.png"), figS3, width = 11, height = 6, dpi = 300, bg = "white")
message("Saved ", file.path(FIG_DIR, "figS3.png"))
