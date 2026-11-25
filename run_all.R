# Run all figure-generation scripts for:
# "Developmental Bias Shapes Evolution on Complex Adaptive Landscapes"
# (Netz, Melo & Milocco)
#
# Run with BLAS threads capped, e.g.:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript run_all.R

scripts <- c(
  "fig2_m_matrices.R",
  "figS3_perlin.R",
  "fig1_linear_maps.R",
  "fig3_basins.R",
  "fig4_fig5_complexity.R",
  "fig6_trajectories.R",
  "motif_example.R"
)

for (s in scripts) {
  message("=== Running ", s, " ===")
  status <- system2("Rscript", s)
  if (status != 0L) stop("Script failed: ", s)
}
message("All figure scripts completed.")
