


### Libraries
library(tidyverse)
library(deSolve)
library(Rlab)
library(gtools)
library(ambient)
library(future)
library(future.apply)
library(rootSolve)
library(viridis)
library(ggnewscale)
library(MASS)

###
Gmat1 <- 5*matrix(c(0.2, 0, 0, 0.2), nrow = 2)
Gmat2 <- 5*matrix(c(0.33333,  0.26666,  0.26666, 0.33333), nrow = 2)
Gmat3 <- 5*matrix(c(0.33333, -0.26666, -0.26666, 0.33333), nrow = 2)

G_matrices <- array(c(Gmat1, Gmat2, Gmat3), dim = c(2, 2, 3))


count_local_maxima <- function(z, x, y) {
  nrow_z <- nrow(z)
  ncol_z <- ncol(z)
  
  z_pad <- matrix(-Inf, nrow = nrow_z + 2, ncol = ncol_z + 2)
  z_pad[2:(nrow_z + 1), 2:(ncol_z + 1)] <- z
  
  maxima <- matrix(FALSE, nrow = nrow_z, ncol = ncol_z)
  
  for (i in 2:(nrow_z + 1)) {
    for (j in 2:(ncol_z + 1)) {
      neighborhood <- z_pad[(i - 1):(i + 1), (j - 1):(j + 1)]
      maxima[i - 1, j - 1] <- z_pad[i, j] == max(neighborhood)
    }
  }
  
  idx <- which(maxima, arr.ind = TRUE)
  
  data.frame(
    row = idx[, 1],
    col = idx[, 2],
    x   = x[idx[, 1]],
    y   = y[idx[, 2]],
    value = z[idx]
  )
  
}


### Parameters
pop_size = 1000
G <- 10000
mut_rate <- 0.05 #0.01
mut_dist <- 0.05
threshold <- 500



# Simulation function

run_rep <- function(rep_id, m, x1, x2, my_freq, my_seed, maxima_) {
  #print(rep_id)
  
  genotype <- rep(c(x1,x2) %*% solve(t(G_matrices[, , m])), each=2)/2
  
  pop <- matrix(rep(genotype, pop_size), nrow = pop_size, byrow = TRUE)
  phenotypes <- matrix(rep(c(x1, x2), pop_size), nrow = pop_size, byrow = TRUE)
  
  fitness <-   rep((gen_perlin(
    x = x1,
    y = x2,
    frequency = my_freq,
    seed = my_seed
  ) + 1) / 2, pop_size)
  
  
  results <- data.frame(
    motif  = numeric(),
    rep     = numeric(),
    x1 = numeric(),
    x2 = numeric(),
    peak_reached = numeric(),
    time    = numeric(),
    mean_X1 = numeric(),
    mean_X2 = numeric()
    
  )
  
  
  # generation loop
  for (g in 1:G) {
    
    #sample parents in a weighted lottery
    m_inds <- sample(1:(pop_size/2),
                     size = pop_size,
                     replace = TRUE,
                     prob = (fitness[1:(pop_size/2)] - min(fitness) + 0.001))
    
    f_inds <- sample((pop_size/2):pop_size,
                     size = pop_size,
                     replace = TRUE,
                     prob = (fitness[(pop_size/2):pop_size] - min(fitness) + 0.001))
    
    
    parents <- cbind(m_inds, f_inds)
    
    
    pop <- cbind(
      pop[cbind(parents[,1], sample(1:1, nrow(pop), TRUE))],
      pop[cbind(parents[,2], sample(1:1, nrow(pop), TRUE))],#
      pop[cbind(parents[,1], sample(3:4, nrow(pop), TRUE))],
      pop[cbind(parents[,2], sample(3:4, nrow(pop), TRUE))]#
    )
    
    # Mutation
    mutations <- rbern(pop_size * 4, mut_rate/2)
    mut_indices <- which(as.logical(mutations)) - 1
    
    ind <- floor(mut_indices / 4) + 1
    loc <- (mut_indices %% 4) + 1
    pop[cbind(ind, loc)] <- pop[cbind(ind, loc)] + rnorm(length(ind), mean = 0, sd = mut_dist)
    
    
    P <- cbind(
      pop[,1] + pop[,2],
      pop[,3] + pop[,4]
    )
    
    phenotypes <- P %*% t(G_matrices[, , m])
    
    trait1 <- phenotypes[,1]
    trait2 <- phenotypes[,2]
    
    fitness <- ((gen_perlin(trait1, trait2, frequency = my_freq, seed = my_seed) + 1) / 2)
    
    
    if (g %% 10 == 0) {
      df <- as.data.frame( phenotypes)
      names(df) <- c("x1", "x2")
      
      condition <-   sapply(seq_len(nrow(maxima_)), function(i) {
        sum(
          dplyr::between(df$x1, maxima_$xmin[i], maxima_$xmax[i]) &
            dplyr::between(df$x2, maxima_$ymin[i], maxima_$ymax[i])
        ) >= threshold
      })
      
      if (any(condition)) {
        results <- bind_rows(
          results,
          tibble(
            motif  = m,
            rep     = rep_id,
            x1 = x1,
            x2 = x2,
            peak_reached = which(condition)[1],
            time    = g,
            mean_X1 = mean(df$x1),
            mean_X2 = mean(df$x2)
            
          )
        )
        
        break
      }
    }
    if (g == G) {
      df <- as.data.frame( phenotypes)
      names(df) <- c("x1", "x2")
      
      results <- bind_rows(
        results,
        tibble(
          motif  = m,
          rep     = rep_id,
          x1 = x1,
          x2 = x2,
          peak_reached = 0,
          time    = g,
          mean_X1 = mean(df$x1),
          mean_X2 = mean(df$x2)
          
        )
      )
    }
  }
  results
}



# Simulation loop
all_results <- list()
k <- 1

targets1 <- rep(seq(2, 8, length.out = 31), each = 31)
targets2 <- rep(seq(0, 6, length.out = 31), 31)

my_freq <- 0.15
my_seed <- 678


res <- 500 # Resolution (100x100)

x_coords <- seq(-1, 15, length.out = res)
y_coords <- seq(-1, 15, length.out = res)
noise_matrix <- outer(x_coords, y_coords, function(x, y) {
  gen_perlin(x, y, frequency = my_freq, seed = my_seed)
})

maxima <- count_local_maxima(noise_matrix, x_coords, y_coords)

maxima <- dplyr::filter(maxima, row > 1, row < 500, col > 1, col < 500)


maxima_ <- data.frame(
  xmin = maxima$x - 0.1,
  xmax = maxima$x + 0.1,
  ymin = maxima$y - 0.1,
  ymax = maxima$y + 0.1
)


plan(multisession)   # works on Windows/Mac/Linux
options(future.rng.onMisuse = "error")

for (m in 1:3) {
  print(paste0("motif: ", m))
  
  for (repl_ in 1:5) {
    print(paste0("repl: ", repl_))
    
    #for(i_ in 1:length(targets1)){
    #  run_rep(repl_, m, targets1[i_], targets2[i_], my_freq, my_seed, maxima_)
    #}
    
    df_list <- future_lapply(seq_along(targets1), function(i_)
      run_rep(repl_, m, targets1[i_], targets2[i_], my_freq, my_seed, maxima_), future.seed = TRUE)
    results_ <- do.call(rbind, df_list)
    
    all_results[[k]] <- results_
    k <- k + 1
  }
}

Sys.time()


final_results <- do.call(rbind, all_results)

# Convert matrix to data frame
noise_df <- expand.grid(x = x_coords, y = y_coords)
noise_df$z <- as.vector(noise_matrix)

# Add columns for the closest maxima
final_results$x1_final <- NA
final_results$x2_final <- NA
final_results$peak_reached2 <- NA

for (i in seq_len(nrow(final_results))) {
  # Calculate Euclidean distances to all maxima
  distances <- sqrt((maxima$x - final_results$mean_X1[i])^2 +
                      (maxima$y - final_results$mean_X2[i])^2)
  
  # Find the index of the closest maxima
  closest <- which.min(distances)
  
  # Assign the corresponding maxima values
  final_results$x1_final[i] <- maxima$x[closest]
  final_results$x2_final[i] <- maxima$y[closest]
  final_results$peak_reached2[i] <- closest
}

# which haven't reached?
filter(final_results, time == 10000)
# distances from optimum assigned
final_results %>% mutate(
  distance = sqrt((mean_X1 - x1_final)^2 + (mean_X2 - x2_final)^2)
) %>%
  arrange(desc(distance))


# determine frequency of most common outcome
processed_results <- final_results %>%
  dplyr::count(motif, x1, x2, peak_reached2) %>%
  group_by(motif, x1, x2) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

processed_results$n <- processed_results$n / 10


### SVM border inference
data <- dplyr::select(final_results, motif, rep, x1, x2, peak_reached2)


library(e1071)
library(dplyr)

extract_pixel_borders <- function(landscape_data, res = 200) {
  # 1. Fit SVM
  landscape_data$peak_reached <- as.factor(landscape_data$peak_reached)
  svm_model <- svm(
    peak_reached ~ x1 + x2,
    data = landscape_data,
    kernel = "radial",
    cost = 10,
    gamma = 0.5,
    cross = 10
  )
  
  print(svm_model$tot.accuracy)
  
  # 2. Create high-res grid
  x1_range <- seq(min(landscape_data$x1),
                  max(landscape_data$x1),
                  length.out = res)
  x2_range <- seq(min(landscape_data$x2),
                  max(landscape_data$x2),
                  length.out = res)
  grid_points <- expand.grid(x1 = x1_range, x2 = x2_range)
  
  # 3. Predict and reshape to matrix
  grid_points$peak_num <- as.numeric(predict(svm_model, grid_points))
  grid_mat <- matrix(grid_points$peak_num, nrow = res, ncol = res)
  
  # 4. Shift-and-Compare Logic (Edge Detection)
  # Check horizontal differences
  diff_h <- grid_mat[-1, ] != grid_mat[-res, ]
  # Check vertical differences
  diff_v <- grid_mat[, -1] != grid_mat[, -res]
  
  # Combine into a logical matrix
  borders <- matrix(FALSE, nrow = res, ncol = res)
  borders[-res, ] <- borders[-res, ] | diff_h
  borders[, -res] <- borders[, -res] | diff_v
  
  # 5. Extract coordinates of the 'TRUE' border pixels
  border_indices <- which(borders, arr.ind = TRUE)
  
  if (nrow(border_indices) == 0)
    return(data.frame())
  
  data.frame(x1 = x1_range[border_indices[, 1]], x2 = x2_range[border_indices[, 2]])
}


# List to store results
all_borders_list <- list()
landscape_ids <- unique(data$motif) # Ensure this column exists

for (id in landscape_ids) {
  message("Processing Landscape: ", id)
  
  # Subset and run function
  sub_data <- data[data$motif == id, ]
  borders <- extract_pixel_borders(sub_data, res = 200)
  
  if (nrow(borders) > 0) {
    borders$motif <- id
    all_borders_list[[as.character(id)]] <- borders
  }
}

# Combine into the final dataframe
master_border_df <- do.call(rbind, all_borders_list)



selected_motif = 3

plot_motif3 <- ggplot(noise_df, aes(x = x, y = y)) +
  
  # --- Background: Perlin noise (continuous fill)
  geom_raster(aes(fill = (z + 1) / 2)) +
  scale_fill_viridis(
    name = "Fitness",
    option = "viridis",
    direction = 1,
    #guide = "none"
  ) +
  
  # --- Overlay: domain borders only (discrete fill)
  geom_segment(
    data = dplyr::filter(final_results, motif == selected_motif),
    aes(
      x = x1,
      y = x2,
      xend = x1_final,
      yend = x2_final
    ),
    arrow = arrow(length = unit(0.05, "cm")),
    color = "orange",
    linewidth = 0.1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = dplyr::filter(final_results, motif == selected_motif),
    aes(x = x1_final, y = x2_final),
    color = "red",
    size = 1,
    shape = 4,
    inherit.aes = FALSE
  ) +
  new_scale_fill() +
  geom_point(
    data = dplyr::filter(processed_results, motif == selected_motif),
    aes(x = x1, y = x2, color = n),
    size = 0.2,
    inherit.aes = FALSE
  ) +
  scale_colour_gradient(low = "white", high = "black") +
  geom_tile(data = dplyr::filter(master_border_df, motif == selected_motif), aes(x = x1, y = x2), fill = "white") +
  geom_tile(data = dplyr::filter(master_border_df, motif == 1), aes(x = x1, y = x2), fill = "red") +

  coord_equal(xlim = c(0, 10),
              ylim = c(0, 7),
              clip = "on") +
  
  #facet_wrap(~ motif,ncol = 2, nrow = 2) +
  theme_minimal() +
  theme(strip.text = element_blank()) +
  labs(
    title = "",
    color = "Max freq.",
    x = expression("Phenotype" ~ x[1]),
    y = expression("Phenotype" ~ x[2])
  )





####
# Evolutionary simulations with a single peak, 3 G-matrices, time series


### Parameters
pop_size = 1000
G <- 5000
mut_rate <- 0.05
mut_dist <- 0.05



# Simulation function

run_rep <- function(rep_id, m, x1, x2, my_freq, my_seed) {
  #print(rep_id)
  
  genotype <- rep(c(x1,x2) %*% solve(t(G_matrices[, , m])), each=2)/2
  
  pop <- matrix(rep(genotype, pop_size), nrow = pop_size, byrow = TRUE)
  phenotypes <- matrix(rep(c(x1, x2), pop_size), nrow = pop_size, byrow = TRUE)

  
  fitness <- rep(exp(-((x1 -5.5)^2 + (x2 -5.5)^2)/(2*3^2)),pop_size)
  
  results <- data.frame(
    time    = numeric(0),
    mean_x1 = numeric(0),
    mean_x2 = numeric(0),
    Gvar1   = numeric(0),
    Gvar2   = numeric(0),
    Gcov = numeric(0)
  )
  
  
  # generation loop
  for (g in 0:G) {
    
    #sample parents in a weighted lottery
    m_inds <- sample(1:(pop_size/2),
                     size = pop_size,
                     replace = TRUE,
                     prob = (fitness[1:(pop_size/2)] - min(fitness) + 0.001))
    
    f_inds <- sample((pop_size/2):pop_size,
                     size = pop_size,
                     replace = TRUE,
                     prob = (fitness[(pop_size/2):pop_size] - min(fitness) + 0.001))
    
    
    parents <- cbind(m_inds, f_inds)
    
    
    pop <- cbind(
      pop[cbind(parents[,1], sample(1:1, nrow(pop), TRUE))],
      pop[cbind(parents[,2], sample(1:1, nrow(pop), TRUE))],#
      pop[cbind(parents[,1], sample(3:4, nrow(pop), TRUE))],
      pop[cbind(parents[,2], sample(3:4, nrow(pop), TRUE))]#
    )
    
    # Mutation
    mutations <- rbern(pop_size * 4, mut_rate/2)
    mut_indices <- which(as.logical(mutations)) - 1
    
    ind <- floor(mut_indices / 4) + 1
    loc <- (mut_indices %% 4) + 1
    pop[cbind(ind, loc)] <- pop[cbind(ind, loc)] + rnorm(length(ind), mean = 0, sd = mut_dist)
    

    P <- cbind(
      pop[,1] + pop[,2],
      pop[,3] + pop[,4]
    )
    
    phenotypes <- P %*% t(G_matrices[, , m])
    
    trait1 <- phenotypes[,1]
    trait2 <- phenotypes[,2]
    
    #fitness <- ((gen_perlin(trait1, trait2, frequency = my_freq, seed = my_seed) + 1) / 2)
    
    fitness <- exp(-((trait1 -5.5)^2 + (trait2 -5.5)^2)/(2*3^2))
    
    if (g %% 10 == 0) {
      output_data <- as.data.frame( phenotypes)
      names(output_data) <- c("x1", "x2")
      
      
      results <- rbind(
        results,
        data.frame(
          time    = g,
          mean_x1 = mean(output_data$x1),
          mean_x2 = mean(output_data$x2),
          Gvar1   = var(output_data$x1),
          Gvar2   = var(output_data$x2),
          Gcov    = cov(output_data$x1, output_data$x2)
        )
      )
    }
    
  }
  results$motif <- m
  results$rep_id <- rep_id
  results$origin_x <- x1
  results$origin_y <- x2
  
  
  results
}


all_results <- list()
k <- 1

targets1 <- c(1,10,5.5,10)
targets2 <- c(5.5,1,10,10)

my_freq <- 0.1
my_seed <- 343
res <- 500 # Resolution 

x_coords <- seq(0, 15, length.out = res)
y_coords <- seq(0, 15, length.out = res)
#noise_matrix <- outer(x_coords, y_coords, function(x, y) {
#  gen_perlin(x, y, frequency = my_freq, seed = my_seed)
#})

noise_matrix <- outer(x_coords, y_coords, function(x, y) {
  exp(-((x -5.5)^2 + (y -5.5)^2)/(2*3^2))
})

maxima <- count_local_maxima(noise_matrix, x_coords, y_coords)

maxima <- dplyr::filter(maxima, row > 1, row < 500, col > 1, col < 500)


plan(multisession)   # works on Windows/Mac/Linux
options(future.rng.onMisuse = "error")

for (m in 1:3) {
  print(paste0("motif: ", m))
  
  for (repl_ in 1:1) {
    print(paste0("repl: ", repl_))
    
    #for(i_ in 1:length(targets1)){
    #  run_rep(repl_, m, targets1[i_], targets2[i_], my_freq, my_seed, maxima_)
    #}
    
    df_list <- future_lapply(seq_along(targets1), function(i_)
      run_rep(repl_, m, targets1[i_], targets2[i_], my_freq, my_seed), future.seed = TRUE)
    results_ <- do.call(rbind, df_list)
    
    all_results[[k]] <- results_
    k <- k + 1
  }
}

Sys.time()


trajectories <- do.call(rbind, all_results)

noise_df3 <- expand.grid(x = x_coords, y = y_coords)
noise_df3$z <- as.vector(noise_matrix)

plot1 <- ggplot(noise_df3, aes(x = x, y = y)) +
  
  # --- Background: Perlin noise (continuous fill)
  geom_raster(aes(fill = z)) +
  scale_fill_viridis(
    name = "Fitness",
    option = "viridis",
    direction = 1,
    #guide = "none"
  ) +
  geom_path(data = trajectories, aes(mean_x1, mean_x2, colour = as.factor(motif), 
                                     group = interaction(motif, origin_x, origin_y)), linewidth = 1.0)+
  geom_point(
    data = maxima,
    aes(x = x, y = y),
    color = "red",
    size = 1,
    shape = 4,  stroke = 1.2,
    inherit.aes = FALSE
  ) +
  scale_colour_discrete(labels = c(" ", " ", " "))+
  coord_equal(xlim = c(0.5, 11),
              ylim = c(0.5, 11),
              clip = "on") +
  labs(
    title = "",
    color = "M-matrix",
    x = expression("Phenotype" ~ x[1]),
    y = expression("Phenotype" ~ x[2])
  )
plot1


plot_motif1 <- plot_motif1 + theme(legend.position = "right")
# Extract legend
legend <- get_legend(
  plot_motif1
)
plot_motif1 <- plot_motif1 + theme(legend.position = "none", plot.margin = margin(t = 25, r = 0, b = 0, l = 0))
plot_motif2 <- plot_motif2 + theme(legend.position = "none",plot.margin = margin(t = 25, r = 0, b = 0, l = 0))
plot_motif3 <- plot_motif3 + theme(legend.position = "none", plot.margin = margin(t = 25, r = 0, b = 0, l = 0))

library(cowplot)

plot1 <- plot1+ guides(fill = "none")

grid <- plot_grid(plot1, plot_motif1, plot_motif2, plot_motif3, ncol = 2, labels = c("A", "B", "C", "D"))

grid

# Add legend on side
final <- plot_grid(
  grid,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.09)
)
final
ggsave("fig1_new.png", height = 8, width = 13, bg="white")

