# Figs 3 and 4



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


###
#Hether & Hohenlohe

motif_1 <- function(a1, a2, g = 1, th = 5) {
  x1 <- a1 / g
  x2 <- a2 * (x1 / (th + x1)) / g
  
  return(c(x1, x2))
}

motif_2 <- function(a1, a2, g = 1, th = 5) {
  x1 <- a1 / g
  x2 <- (a2 * (1 - ((a1 / g) / (th + ((a1) / (g)
  ))))) / (g)
  
  return(c(x1, x2))
}

motif_3 <- function(a1, a2, g = 1, th = 5) {
  x1 <- (2 * th * a1 * sqrt(g)) / ((-a1 + th * g) * sqrt(g) + sqrt(4 * a1 *
                                                                     a2 * g + ((a1 + th * g)^2) * g))
  x2 <- (-(a1 + th * g) * sqrt(g) + sqrt(4 * a1 * a2 * g + ((a1 + th * g)^2) *
                                           g)) / (2 * g * sqrt(g))
  
  return(c(x1, x2))
}

motif_4 <- function(a1, a2, g = 1, th = 5) {
  x1 <- (a1 - a2 - g * th + sqrt(4 * a1 * g * th + (a1 - a2 - g * th)^2)) /
    (2 * g)
  x2 <- (-a1 + a2 - g * th + sqrt(4 * a1 * g * th + (a1 - a2 - g * th)^2)) /
    (2 * g)
  
  return(c(x1, x2))
}

motif_5 <- function(a1, a2, g = 1, th = 5) {
  x1 <- (a1 * a2 - (g^2) * (th^2)) / (g * (a2 + g * th))
  x2 <- (a1 * a2 - (g^2) * (th^2)) / (g * (a1 + g * th))
  
  return(c(x1, x2))
}

motif_6 <- function(a1, a2, g = 1, th = 5) {
  x1 <- a1 / g
  x2 <- a2 / g
  
  return(c(x1, x2))
}

motifsGtoP <- list(
  motif_1  = motif_1,
  motif_2  = motif_2,
  motif_3  = motif_3,
  motif_4  = motif_4,
  motif_5  = motif_5,
  motif_6  = motif_6
)


motifPG_1 <- function(x1, x2, g = 1, th = 5) {
  a1 <- g * x1
  a2 <- (g * x2 * (th + x1)) / (x1)
  return(c(a1, a2))
}

motifPG_2 <- function(x1, x2, g = 1, th = 5) {
  a1 <- g * x1
  a2 <- (g * x2) / (1 - (x1 / (th + x1)))
  return(c(a1, a2))
}

motifPG_3 <- function(x1, x2, g = 1, th = 5) {
  a1 <- (g * x1) / (1 - ((x2) / (th + x2)))
  a2 <- (g * x2 * (th + x1)) / (x1)
  
  return(c(a1, a2))
}

motifPG_4 <- function(x1, x2, g = 1, th = 5) {
  a1 <- (g * x1 * (th + x2)) / (th)
  a2 <- (g * x2 * (th + x1)) / (th)
  return(c(a1, a2))
}

motifPG_5 <- function(x1, x2, g = 1, th = 5) {
  a1 <- (g * x1 * (th + x2)) / (x2)
  a2 <- (g * x2 * (th + x1)) / (x1)
  return(c(a1, a2))
}

motifPG_6 <- function(x1, x2, g = 1, th = 5) {
  a1 <- x1 * g
  a2 <- x2 * g
  return(c(a1, a2))
}

motifsPtoG <- list(
  motifPG_1  = motifPG_1,
  motifPG_2  = motifPG_2,
  motifPG_3  = motifPG_3,
  motifPG_4  = motifPG_4,
  motifPG_5  = motifPG_5,
  motifPG_6  = motifPG_6
)


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
mut_rate <- 0.01
mut_dist <- 0.05
threshold <- 500

# Simulation function

run_rep <- function(rep_id, m, x1, x2, my_freq, my_seed, maxima_) {
  #print(rep_id)
  
  genotype <- motifsPtoG[[m]](x1, x2)
  
  pop <- replicate(pop_size, genotype, simplify = FALSE)
  phenotypes <- replicate(pop_size, c(x1, x2), simplify = FALSE)
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
    #sample offspring in a weighted lottery
    inds <- sample(
      1:pop_size,
      size = pop_size,
      replace = TRUE,
      prob = (fitness - min(fitness) + 0.001)
    )
    
    pop <- pop[inds]
    phenotypes <- phenotypes[inds]
    fitness <- fitness[inds]
    
    # Mutation
    mutations <- rbern(pop_size * 2, mut_rate)
    mut_indices <- which(as.logical(mutations)) - 1
    
    for (mut in mut_indices) {
      ind <- floor(mut / 2) + 1
      loc <- (mut %% 2) + 1
      
      pop[[ind]][loc] <- max(0.0, pop[[ind]][loc] + rnorm(1, 0, mut_dist))
      
      phenotypes[[ind]] <- motifsGtoP[[m]](pop[[ind]][[1]], pop[[ind]][[2]])
      
      #Updating fitness
      fitness[[ind]] <- (
        gen_perlin(
          x = phenotypes[[ind]][1],
          y = phenotypes[[ind]][2],
          frequency = my_freq,
          seed = my_seed
        ) + 1
      ) / 2
      
    }
    
    
    
    if (g %% 10 == 0) {
      df <- as.data.frame(do.call(rbind, phenotypes))
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
            peak_reached = which(condition),
            time    = g,
            mean_X1 = mean(df$x1),
            mean_X2 = mean(df$x2)
            
          )
        )
        
        break
      }
    }
    if (g == G) {
      df <- as.data.frame(do.call(rbind, phenotypes))
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

targets1 <- rep(seq(2, 10, length.out = 41), each = 41)
targets2 <- rep(seq(2, 10, length.out = 41), 41)

#my_seed <- 1234
my_seed <- 1245
#my_freq <- 0.1
my_freq <- 0.5
res <- 500 # Resolution (100x100)

x_coords <- seq(0, 15, length.out = res)
y_coords <- seq(0, 15, length.out = res)
noise_matrix <- outer(x_coords, y_coords, function(x, y) {
  gen_perlin(x, y, frequency = my_freq, seed = my_seed)
})

maxima <- count_local_maxima(noise_matrix, x_coords, y_coords)

maxima <- filter(maxima, row > 1, row < 500, col > 1, col < 500)


maxima_ <- data.frame(
  xmin = maxima$x - 0.1,
  xmax = maxima$x + 0.1,
  ymin = maxima$y - 0.1,
  ymax = maxima$y + 0.1
)


plan(multisession)   # works on Windows/Mac/Linux
options(future.rng.onMisuse = "error")

for (m in 1:6) {
  print(paste0("motif: ", m))
  
  for (repl_ in 1:10) {
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
final_results2$x1_final <- NA
final_results2$x2_final <- NA
final_results2$peak_reached2 <- NA

for (i in seq_len(nrow(final_results2))) {
  # Calculate Euclidean distances to all maxima
  distances <- sqrt((maxima$x - final_results2$mean_X1[i])^2 +
                      (maxima$y - final_results2$mean_X2[i])^2)
  
  # Find the index of the closest maxima
  closest <- which.min(distances)
  
  # Assign the corresponding maxima values
  final_results2$x1_final[i] <- maxima$x[closest]
  final_results2$x2_final[i] <- maxima$y[closest]
  final_results2$peak_reached2[i] <- closest
}

# which haven't reached?
filter(final_results2, time == 10000)
# distances from optimum assigned
final_results2 %>% mutate(
  distance = sqrt((mean_X1 - x1_final)^2 + (mean_X2 - x2_final)^2)
) %>%
  arrange(desc(distance))


# determine frequency of most common outcome
processed_results2 <- final_results2 %>%
  dplyr::count(motif, x1, x2, peak_reached2) %>%
  group_by(motif, x1, x2) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

processed_results2$n <- processed_results2$n / 10


### SVM border inference
data <- dplyr::select(final_results2, motif, rep, x1, x2, peak_reached2)


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
    gamma = 0.5
  )
  
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


ggplot(noise_df, aes(x = x, y = y)) +
  
  # --- Background: Perlin noise (continuous fill)
  geom_raster(aes(fill = (z + 1) / 2)) +
  scale_fill_viridis(
    name = "fitness",
    option = "viridis",
    direction = 1,
    #guide = "none"
  ) +
  
  # --- Overlay: domain borders only (discrete fill)
  geom_segment(
    data = final_results2,
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
    data = final_results2,
    aes(x = x1_final, y = x2_final),
    size = 1,
    shape = 4,
    color = "red",
    inherit.aes = FALSE
  ) +
  geom_tile(data = select(filter(master_border_df, motif == 6), -motif),
            aes(x = x1, y = x2),
            fill = "red") +
  new_scale_fill() +
  geom_point(
    data = processed_results2,
    aes(x = x1, y = x2, color = n),
    size = 0.2,
    inherit.aes = FALSE
  ) +
  scale_colour_gradient(low = "white", high = "black") +
  geom_tile(data = master_border_df, aes(x = x1, y = x2), fill = "white") +
  
  coord_equal(xlim = c(0, 11.5),
              ylim = c(0, 11.5),
              clip = "on") +
  
  facet_wrap( ~ motif, labeller = labeller(
    motif = function(x)
      paste("motif", x)
  )) +
  theme_minimal() +
  labs(
    title = "",
    color = "max peak frequency",
    x = expression(x[1]),
    y = expression(x[2])
  )


ggsave("fig5.png", width = 10, height = 7)


###
# Supplementary figure: trace out analytical separatrices

library(ambient)
library(ggplot2)
library(viridis)
library(future.apply)


determine_manifolds <- function(j) {
  Sigma <- matrix(c(5, 0, 0, 5), nrow = 2)
  J_fun <- J_matrix[[j]]

  
  G_fun <- function(x, y) {
    J_fun(x, y, th = 5) %*% Sigma %*% t(J_fun(x, y, th = 5))
  }
  

  # Bilinear Interpolation
  bilinear_interp <- function(x, y, mat) {
    # Normalize to grid indices
    ix <- (x - domain_min) / (domain_max - domain_min) * (length(x_coords) - 1) + 1
    iy <- (y - domain_min) / (domain_max - domain_min) * (length(y_coords) - 1) + 1
    
    # Clamp indices to be within the valid grid range
    i1 <- pmax(pmin(floor(ix), res - 1), 1)
    j1 <- pmax(pmin(floor(iy), res - 1), 1)
    
    # Prevent i2 and j2 from exceeding res by checking if i1 and j1 are at max index
    i2 <- pmin(i1 + 1, res - 1)
    j2 <- pmin(j1 + 1, res - 1)
    
    # Calculate interpolation weights
    wx <- ix - i1
    wy <- iy - j1
    
    # Safely index the matrix using the clamped indices
    #print(c(i1, j1, i2, j2))
    f11 <- mat[cbind(i1, j1)]
    f21 <- mat[cbind(i2, j1)]
    f12 <- mat[cbind(i1, j2)]
    f22 <- mat[cbind(i2, j2)]
    
    # Return interpolated value
    return((1 - wx) * (1 - wy) * f11 +
             wx * (1 - wy) * f21 +
             (1 - wx) * wy * f12 +
             wx * wy * f22)
  }
  interp_grad <- function(x, y) {
    c(bilinear_interp(x, y, grad_x),
      bilinear_interp(x, y, grad_y))
  }
  
  interp_hessian <- function(x, y) {
    hxx <- bilinear_interp(x, y, dxx)
    hyy <- bilinear_interp(x, y, dyy)
    hxy <- bilinear_interp(x, y, dxy)
    matrix(c(hxx, hxy, hxy, hyy), 2)
  }
  

  # Newton Critical Point Finder
  find_critical <- function(x0,
                            y0,
                            max_iter = 50,
                            tol = 1e-6) {
    x <- x0
    y <- y0
    
    for (i in 1:max_iter) {
      g <- interp_grad(x, y)
      if (sqrt(sum(g^2)) < tol)
        break
      
      H <- interp_hessian(x, y)
      det_H <- det(H)
      
      if (abs(det_H) < 1e-6) {
        # Skip poorly conditioned Hessians (potential saddle failures)
        return(NULL)
      }
      
      step <- tryCatch(
        solve(H, g),
        error = function(e)
          return(NULL)
      )
      if (is.null(step))
        return(NULL)
      
      # Update the positions
      x <- x - step[1]
      y <- y - step[2]
      
      # Prevent going out of bounds
      if (x < domain_min ||
          x > domain_max || y < domain_min || y > domain_max)
        return(NULL)
    }
    
    return(c(x, y))
  }
  # Seed random initial guesses
  set.seed(1)
  seed_spacing <- 0.5   # adjust
  seed_x <- seq(domain_min, domain_max, by = seed_spacing)
  seed_y <- seq(domain_min, domain_max, by = seed_spacing)
  
  seed_grid <- expand.grid(seed_x, seed_y)
  
  candidates <- lapply(1:nrow(seed_grid), function(i) {
    find_critical(seed_grid[i, 1], seed_grid[i, 2])
  })
  
  crit_pts <- unique(do.call(rbind, candidates))
  crit_pts <- crit_pts[complete.cases(crit_pts), ]
  

  # Classify critical points
  classify_point <- function(pt) {
    H <- interp_hessian(pt[1], pt[2])
    H <- H + 1e-8 * diag(2)
    
    if (any(!is.finite(H)))
      return(NA)
    
    Gmat <- G_fun(pt[1], pt[2])
    J <- Gmat %*% H
    
    if (any(!is.finite(J)))
      return(NA)
    
    eig <- tryCatch(
      eigen(J),
      error = function(e)
        return(NA)
    )
    if (is.na(eig)[1])
      return(NA)
    
    sum(Re(eig$values) > 0)
  }
  
  types <- apply(crit_pts, 1, classify_point)
  saddles <- crit_pts[types == 1, , drop = FALSE]
  
  cat("Number of saddles detected:", nrow(saddles), "\n")
  
  # Flow Field + RK4 Integrator
  flow_field_unit <- function(z, eps = 1e-10) {
    g <- interp_grad(z[1], z[2])
    Gmat <- G_fun(z[1], z[2])
    v <- as.numeric(Gmat %*% g)
    
    norm_v <- sqrt(sum(v^2))
    v / (norm_v + eps)
  }
  
  rk4_step_unit <- function(z, ds) {
    k1 <- flow_field_unit(z)
    k2 <- flow_field_unit(z + 0.5 * ds * k1)
    k3 <- flow_field_unit(z + 0.5 * ds * k2)
    k4 <- flow_field_unit(z + ds * k3)
    
    z + ds * (k1 + 2 * k2 + 2 * k3 + k4) / 6
  }
  
  
  trace_branch <- function(z0,
                           dir,
                           ds = 0.01,
                           max_steps = 4000) {
    z <- z0 + 1e-6 * dir
    path <- matrix(NA, max_steps, 2)
    
    for (i in 1:max_steps) {
      path[i, ] <- z
      z <- rk4_step_unit(z, ds)
      
      if (any(z < domain_min | z > domain_max))
        break
    }
    
    na.omit(path)
  }
  
  
  # Compute stable manifolds and store in dataframe
  manifold_df_list <- list()
  row_counter <- 1
  
  for (i in 1:nrow(saddles)) {
    z <- saddles[i, ]
    H <- interp_hessian(z[1], z[2])
    Gmat <- G_fun(z[1], z[2])
    J <- Gmat %*% H
    
    
    if (any(!is.finite(J)))
      next
    
    eig <- tryCatch(
      eigen(J),
      error = function(e)
        NULL
    )
    
    if (is.null(eig))
      next
    
    #eig <- eigen(J)
    vals <- Re(eig$values)
    vecs <- Re(eig$vectors)
    
    stable_idx <- which(vals < 0)
    if (length(stable_idx) != 1)
      next
    
    stable_vec <- vecs[, stable_idx]
    stable_vec <- stable_vec / sqrt(sum(stable_vec^2))
    
    branch1 <- trace_branch(z, stable_vec, ds = -0.01)
    branch2 <- trace_branch(z, -stable_vec, ds = -0.01)
    
    for (b in list(branch1, branch2)) {
      if (is.null(b) || nrow(b) == 0)
        next
      
      df_branch <- data.frame(
        x = b[, 1],
        y = b[, 2],
        motif = j,
        saddle_id = i,
        branch_id = row_counter,
        step = seq_len(nrow(b))
      )
      
      manifold_df_list[[length(manifold_df_list) + 1]] <- df_branch
      row_counter <- row_counter + 1
    }
  }
  motif_df <- do.call(rbind, manifold_df_list)
  motif_df
}


all_manifolds <- list()
manifolds_df <- data.frame()


domain_min <- 0
domain_max <- 12

dx_target <- 0.01
res <- ceiling((domain_max - domain_min) / dx_target)

x_coords <- seq(0, 12, length.out = res)
y_coords <- seq(0, 12, length.out = res)

dx <- x_coords[2] - x_coords[1]
dy <- y_coords[2] - y_coords[1]

# Perlin landscape
# simple
# noise_matrix <- outer(x_coords, y_coords, function(x, y) {
#    gen_perlin(x, y, frequency = 0.1, seed = 1234)
#  })

# complex
noise_matrix <- outer(x_coords, y_coords, function(x, y) {
  gen_perlin(x, y, frequency = 0.5, seed = 1245)
})


# Gradient (central difference)
grad_x <- rbind(
  (noise_matrix[2, ] - noise_matrix[1, ]) / dx,
  (noise_matrix[3:res, ] - noise_matrix[1:(res - 2), ]) / (2 * dx),
  (noise_matrix[res, ] - noise_matrix[res - 1, ]) / dx
)

grad_y <- cbind(
  (noise_matrix[, 2] - noise_matrix[, 1]) / dy,
  (noise_matrix[, 3:res] - noise_matrix[, 1:(res - 2)]) / (2 * dy),
  (noise_matrix[, res] - noise_matrix[, res - 1]) / dy
)

# Hessian components
dxx <- rbind((grad_x[2, ] - grad_x[1, ]) / dx,
             (grad_x[3:res, ] - grad_x[1:(res - 2), ]) / (2 * dx),
             (grad_x[res, ] - grad_x[res - 1, ]) / dx)

dyy <- cbind((grad_y[, 2] - grad_y[, 1]) / dy,
             (grad_y[, 3:res] - grad_y[, 1:(res - 2)]) / (2 * dy),
             (grad_y[, res] - grad_y[, res - 1]) / dy)

dxy <- (rbind(grad_y[3:res, ], grad_y[res, ]) -
          rbind(grad_y[1, ], grad_y[1:(res - 2), ])) / (2 * dx)


plan(multisession)

all_manifolds <- future_lapply(1:6, function(j)
  determine_manifolds(j), future.seed = TRUE)

manifolds_df <- do.call(rbind, all_manifolds)

restricted_df <- subset(manifolds_df, x >= 2 & x <= 10 &
                          y >= 2 & y <= 10)

ggplot(noise_df, aes(x = x, y = y)) +
  
  # --- Background: Perlin noise (continuous fill)
  geom_raster(aes(fill = (z + 1) / 2)) +
  scale_fill_viridis(
    name = "fitness",
    option = "viridis",
    # try "magma", "inferno", "plasma", "cividis"
    direction = 1,
    #guide = "none"
  ) +
  
  # --- Reset fill scale
  #new_scale_fill() +
  
  # --- Overlay: domain borders only (discrete fill)
  geom_segment(
    data = final_results2,
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
    data = final_results2,
    aes(x = x1_final, y = x2_final),
    size = 1,
    shape = 4,
    color = "red",
    inherit.aes = FALSE
  ) +
  #geom_tile(data = select(filter(master_border_df, motif == 6), -motif), aes(x = x1, y = x2), fill = "red")+
  
  new_scale_fill() +
  geom_point(
    data = processed_results2,
    aes(x = x1, y = x2, color = n),
    size = 0.2,
    inherit.aes = FALSE
  ) +
  scale_colour_gradient(low = "white", high = "black") +
  geom_tile(data = master_border_df, aes(x = x1, y = x2), fill = "white") +
  geom_path(
    data = restricted_df,
    aes(x, y, group = branch_id),
    color = "black",
    linewidth = 0.2
  ) +
  coord_equal(xlim = c(0, 11.5),
              ylim = c(0, 11.5),
              clip = "on") +
  
  facet_wrap( ~ motif, labeller = labeller(
    motif = function(x)
      paste("motif", x)
  )) +
  theme_minimal() +
  labs(
    title = "",
    color = "max peak frequency",
    x = expression(x[1]),
    y = expression(x[2])
  )

ggsave("figSX.png", width = 10, height = 7)


# Determine the error rate of inferred separatrices, to be done



