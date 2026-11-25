# GRN_Landscapes

R scripts reproducing all simulations and figures for *Developmental Bias Shapes Evolution on Complex Adaptive Landscapes* (Netz, Melo & Milocco), powered by the [`grnland`](https://github.com/diogro/grnland) R package.

## Requirements

Install `grnland` and the plotting/analysis packages:

```r
# remotes::install_github("diogro/grnland")
install.packages(c("ambient", "tidyverse", "patchwork", "e1071",
                   "ggnewscale", "igraph", "ggraph"))
```

## Scripts

| Script | Output | Description |
|---|---|---|
| `common.R` | — | Shared simulation defaults, parallel helpers, SVM separatrix fitting, parent-offspring G-matrix estimator |
| `fig1_linear_maps.R` | `figures/fig1.png` | Linear GP maps ($\mathbf{B}_1, \mathbf{B}_2, \mathbf{B}_3$) on single- and two-peak landscapes |
| `fig2_m_matrices.R` | `figures/fig2.png` | Six GRN motifs, diploid $\mathbf{M}$-matrix tensor fields, and sensitivity vectors |
| `fig3_basins.R` | `figures/fig3.png`, `figS1.png`, `figS2.png` | Basins of attraction on simple (freq. 0.1) and complex (freq. 0.5) Perlin landscapes, SVM borders and Lande separatrices |
| `fig4_fig5_complexity.R` | `figures/fig4.png`, `figures/fig5.png`, `figures/fig5_e_d_sigma0.2.png` | Outcome mismatch and $M_{\max}$ alignment across Perlin frequencies (0.3, 0.6, 1.2, 1.8) at $\sigma = 0.2$ (and `SIGMA=0.05` check); alternative alignment measure $e_d = (\mathbf{d}^\top \mathbf{M} \mathbf{d} / \mathbf{d}^\top \mathbf{d}) / (\mathrm{tr}\,\mathbf{M}/2)$, the relative mutational variance along the divergence vector, with $\mathbf{M}$ at the end or starting point |
| `fig6_trajectories.R` | `figures/fig6.png` | Replicated trajectories of Motifs 1 and 5 ($N = 10{,}000$) tracking parent-offspring $\mathbf{G}$ and local $\mathbf{M}$ |
| `fig6_search.R` | `output/fig6_search_candidates.csv` | Screen for robust Motif 1 vs. 5 divergence from similar initial $\mathbf{M}$-matrices |
| `motif_example.R` | `figures/motif_example.png` | Single-start divergence across all six motifs on a rugged Perlin landscape |
| `figS3_perlin.R` | `figures/figS3.png` | Perlin noise landscapes (frequencies 0.1–1.8) and spatial autocorrelation profiles |
| `trial_selection.R` | `figures/trial_selection_comparison.png` | Sensitivity comparison of `selection_relative_min(0.001)` vs. `selection_proportional()` |
| `run_all.R` | all main & supplementary figures | Master runner script |

## Simulation defaults

All scripts share the diploid sexual population model defined in `common.R`:

- Population size $N = 1{,}000$ diploid individuals ($N = 10{,}000$ in Fig. 6), split evenly between males and females
- Two additive loci ($\alpha_1, \alpha_2$), genotypic value $\alpha_i = \alpha_{i,1} + \alpha_{i,2}$, free recombination
- Per-allele mutation rate $\mu = 0.025$, Gaussian mutational effect SD $\sigma = 0.05$ ($\sigma = 0.2$ in Figs. 4–5), alleles truncated at $0$ for GRNs (unbounded for linear maps in Fig. 1)
- Selection weights $w - \min(w) + 0.001$ (`selection_relative_min(0.001)`); set `SELECTION=proportional` to run with Wright-Fisher proportional fitness weights
- Convergence rule: population mean phenotype remains within Euclidean distance $0.1$ of a local peak for two consecutive 10-generation checkpoints (up to 10,000 generations; 3,000 generations in the Fig. 4–5 complexity sweep)

## Running

Cap BLAS threads before running parallel scripts:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript run_all.R
```

Useful environment variables:

- `QUICK=1`: run reduced grids and replicates into `output/quick/` and `figures/quick/` to verify the pipeline.
- `N_CORES=k`: number of parallel worker processes.
- `REDO=1`: recompute cached `.rds` results in `output/`.
- `SIGMA=0.05`: run `fig4_fig5_complexity.R` at $\sigma = 0.05$.
- `SELECTION=proportional`: run any script under proportional fitness selection.
