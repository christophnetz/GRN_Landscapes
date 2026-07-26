###
# Everything to create figure 6/trajectories

### Libraries
library(tidyverse)
library(deSolve)
library(Rlab)
library(gtools)
library(ambient)
library(future)
library(future.apply)
library(rootSolve)
library(ellipse)


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


### Jacobians
mut_dist <- 0.05

Sigma <- mut_rate * matrix(c(2 * mut_dist, 0, 0, 2 * mut_dist), nrow = 2)

J1 <- function(x1, x2, g = 1, th = 5) {
  a1 <- g * x1
  a2 <- (g * x2 * (th + x1)) / (x1)
  
  return(matrix(c(1 / g, a2 * th / (g * (
    a1 + th
  )^2), 0, a1 / (g * (
    th + a1
  ))), nrow = 2))
}

J2 <- function(x1, x2, g = 1, th = 5) {
  a1 <- g * x1
  a2 <- (g * x2) / (1 - (x1 / (th + x1)))
  
  x1 <- a1 / g
  x2 <- (a2 * (1 - ((a1 / g) / (th + ((a1) / (g)
  ))))) / (g)
  
  
  return(matrix(c(1 / g, -a2 * th / (a1 + g * th)^2, 0, (1 - ((a1 / g) /
                                                                (th + ((a1) / (g)))
  )) / (g)), nrow = 2))
}

J3 <- function(x1, x2, g = 1, th = 5) {
  a1 <- (g * x1) / (1 - ((x2) / (th + x2)))
  a2 <- (g * x2 * (th + x1)) / (x1)
  
  return(matrix(c((th * (a1 + 2 * a2 + g * th + sqrt(
    4 * a1 * a2 + (a1 + g * th)^2
  ))) / (2 * (a2 + g * th) * sqrt(4 * a1 * a2 + (a1 + g * th)^2)),
  ((a1 + 2 * a2 + g * th) / sqrt(4 * a1 * a2 + (a1 + g * th)^2) -
     1) / (2 * g),-(4 * a1^2 * th) / (sqrt(4 * a1 * a2 + (a1 + g * th)^2) * (-a1 +
                                                                               g * th + sqrt(
                                                                                 4 * a1 * a2 + (a1 + g * th)^2
                                                                               ))^2),
  a1 / (g * sqrt(4 * a1 * a2 + (a1 + g * th)^2))
  ), nrow = 2))
}

J4 <- function(x1, x2, g = 1, th = 5) {
  a1 <- (g * x1 * (th + x2)) / (th)
  a2 <- (g * x2 * (th + x1)) / (th)
  
  return(matrix(c(
    (1 + (a1 - a2 + g * th) / (sqrt(
      4 * a1 * g * th + (-a1 + a2 + g * th)^2
    ))) / (2 * g),
    (-1 + (a1 - a2 + g * th) / (sqrt(
      4 * a1 * g * th + (-a1 + a2 + g * th)^2
    ))) / (2 * g),
    (-1 + (-a1 + a2 + g * th) / (sqrt(
      4 * a1 * g * th + (-a1 + a2 + g * th)^2
    ))) / (2 * g),
    (1 + (-a1 + a2 + g * th) / (sqrt(
      4 * a1 * g * th + (-a1 + a2 + g * th)^2
    ))) / (2 * g)
  ), nrow = 2))
}

J5 <- function(x1, x2, g = 1, th = 5) {
  a1 <- (g * x1 * (th + x2)) / (x2)
  a2 <- (g * x2 * (th + x1)) / (x1)
  
  return(matrix(c(
    a2 / (g * (a2 + th)), (th * (a2 + th)) / (g * (a1 + th)^2), (th * (a1 + th)) /
      (g * (a2 + th)^2), a1 / (g * (a1 + th))
  ), nrow = 2))
}

J6 <- function(x1, x2, g = 1, th = 5) {
  a1 <- x1 * g
  a2 <- x2 * g
  
  return(matrix(c(1 / g, 0, 0, 1 / g), nrow = 2))
}


J_matrix <- list(
  J1  = J1,
  J2  = J2,
  J3  = J3,
  J4  = J4,
  J5  = J5,
  J6  = J6
)

mat_to_ellipse <- function(M, centre, level = 0.95) {
  e <- as.data.frame(ellipse(M, centre = centre, level = level))
  colnames(e) <- c("ex", "ey")
  e
}


### Parameters
pop_size = 10000
G <- 10000 # 12000
mut_rate <- 0.05 #0.008
mut_dist <- 0.05  #0.04

replicates <- 1
nr_motives <- 2
interval <- 200
samples_per_run <- G / interval
length_df1 <- G * replicates * nr_motives
length_df2 <- samples_per_run * replicates * nr_motives
results <- data.frame(
  motif  = numeric(length_df1),
  rep     = numeric(length_df1),
  time    = numeric(length_df1),
  mean_x1 = numeric(length_df1),
  mean_x2 = numeric(length_df1)
)

results2 <- data.frame(
  motif  = numeric(length_df2),
  rep     = numeric(length_df2),
  time    = numeric(length_df2),
  mean_x1 = numeric(length_df2),
  mean_x2 = numeric(length_df2),
  Gvar1   = numeric(length_df2),
  Gvar2   = numeric(length_df2),
  Gcov    = numeric(length_df2),
  Gvar1_   = numeric(length_df2),
  Gvar2_   = numeric(length_df2),
  Gcov_    = numeric(length_df2)
)


pop <- list()
phenotypes <- list()
fitness <- vector()



### Fitness
# base contributions


fitness_fun <- function(a, b, beta = 100) {
  ea <- exp(beta * a)
  eb <- exp(beta * b)
  
  (ea * a + eb * b) / (ea + eb)
}



t0 <- Sys.time()
for (m in c(1, 5)) {
  # for (m in 1:nr_motives) {
  print(m)
  
  for (rep in 1:replicates) {
    print(rep)
    
    genotype <- rep(motifsPtoG[[m]](1.12, 6.4), each = 2) / 2
    
    pop <- matrix(rep(genotype, pop_size),
                  nrow = pop_size,
                  byrow = TRUE)
    phenotypes <- replicate(pop_size, c(1.12, 6.4), simplify = FALSE)
    fitness <-   rep(fitness_fun(((
      gen_perlin(
        x = 1.12 - 1.0,
        y = 6.4 + 2.0,
        frequency = 0.05,
        seed = 76
      ) + 1
    ) / 2)^2, ((
      gen_perlin(
        x = 1.12 - 2.5,
        y = 6.4 + 1,
        frequency = 0.05,
        seed = 76
      ) + 1
    ) / 2)^2), pop_size)
    
    
    # generation loop
    for (g in 1:G) {
      #sample parents in a weighted lottery
      m_inds <- sample(
        1:(pop_size / 2),
        size = pop_size,
        replace = TRUE,
        prob = (fitness[1:(pop_size / 2)] - min(fitness) + 0.001)
      )
      
      f_inds <- sample((pop_size / 2):pop_size,
                       size = pop_size,
                       replace = TRUE,
                       prob = (fitness[(pop_size / 2):pop_size] - min(fitness) + 0.001)
      )
      
      
      parents <- cbind(m_inds, f_inds)
      
      
      pop <- cbind(pop[cbind(parents[, 1], sample(1:2, nrow(pop), TRUE))], pop[cbind(parents[, 2], sample(1:2, nrow(pop), TRUE))], pop[cbind(parents[, 1], sample(3:4, nrow(pop), TRUE))], pop[cbind(parents[, 2], sample(3:4, nrow(pop), TRUE))])
      
      # Mutation
      mutations <- rbern(pop_size * 4, mut_rate / 2)
      mut_indices <- which(as.logical(mutations)) - 1
      
      ind <- floor(mut_indices / 4) + 1
      loc <- (mut_indices %% 4) + 1
      pop[cbind(ind, loc)] <- pmax(0.0, pop[cbind(ind, loc)] + rnorm(length(ind), mean = 0, sd = mut_dist))
      
      phenotypes <- Map(motifsGtoP[[m]], pop[, 1] + pop[, 2], pop[, 3] + pop[, 4])
      P <- do.call(rbind, phenotypes)
      
      x1 <- P[, 1]
      y1 <- P[, 2]
      
      a <- ((gen_perlin(
        x1 - 1.0, y1 + 2.0, frequency = 0.05, seed = 76
      ) + 1) / 2)^2
      b <- ((gen_perlin(
        x1 - 2.5, y1 + 1, frequency = 0.05, seed = 76
      ) + 1) / 2)^2
      fitness <- fitness_fun(a, b)
      
      output_data <- as.data.frame(do.call(rbind, phenotypes))
      names(output_data) <- c("x1", "x2")
      
      results[(G) * (rep - 1) + (G) * (replicates) * (m - 1) + g, ] <- c(
        m,
        rep,
        g,
        mean(output_data$x1),
        mean(output_data$x2)
      )
      if (g %% interval == 0) {
        #sample parents in a lottery
        m_inds <- sample(
          1:(pop_size / 2),
          size = pop_size,
          replace = TRUE
        )

        f_inds <- sample((pop_size / 2):pop_size,
                         size = pop_size,
                         replace = TRUE
                         )


        parents <- cbind(m_inds, f_inds)

        mean_parent <- (do.call(rbind, phenotypes[m_inds]) +  do.call(rbind, phenotypes[f_inds]))/2

        pop <- cbind(pop[cbind(parents[, 1], sample(1:2, nrow(pop), TRUE))], pop[cbind(parents[, 2], sample(1:2, nrow(pop), TRUE))], pop[cbind(parents[, 1], sample(3:4, nrow(pop), TRUE))], pop[cbind(parents[, 2], sample(3:4, nrow(pop), TRUE))])

        offspring_phenotypes <- do.call(rbind, Map(motifsGtoP[[m]], pop[, 1] + pop[, 2], pop[, 3] + pop[, 4]))

        df <- as.data.frame(cbind(mean_parent, offspring_phenotypes))
        colnames(df) <- c("mp_trait1", "mp_trait2", "off_trait1", "off_trait2")

        fit <- lm(cbind(off_trait1, off_trait2) ~ mp_trait1  + mp_trait2,
                  data = df)

        B <- t(coef(fit)[-1, ])


        P <- cov(mean_parent)
        G_m <- B %*% P

        results2[samples_per_run * (rep - 1) + samples_per_run * (replicates) * (m - 1) + floor(g/interval), ] <- c(
          m,
          rep,
          g,
          mean(output_data$x1),
          mean(output_data$x2),
          var(output_data$x1),
          var(output_data$x2),
          cov(output_data$x1, output_data$x2),
          G_m[1,1],
          G_m[2,2],
          G_m[1,2]
        )
      }
    }
  }
}

Sys.time() - t0

###
# Plot trajectories and G matrices
###
# Plot A
res <- 500
x_coords <- seq(0, 15, length.out = res)
y_coords <- seq(0, 15, length.out = res)


base <- outer(x_coords, y_coords, function(x, y) {
  gen_perlin(x - 1.0, y + 2.0, frequency = 0.05, seed = 76)
})

base <- (base + 1) / 2

terrain_steep <- base^2

base2 <- outer(x_coords, y_coords, function(x, y) {
  gen_perlin(x - 2.5, y + 1.0, frequency = 0.05, seed = 76)
})

base2 <- (base2 + 1) / 2

terrain_steep2 <- base2^2

softmax2 <- function(a, b, beta = 5) {
  ea <- exp(beta * a)
  eb <- exp(beta * b)
  ea / (ea + eb)
}

w <- softmax2(terrain_steep, terrain_steep2, beta = 100)

terrain <- w * terrain_steep + (1 - w) * terrain_steep2


noise_df <- expand.grid(x = x_coords, y = y_coords)

noise_df$z <- as.vector(terrain)


# Recreate matrix from data frame
z_mat <- matrix(noise_df$z, nrow = res, ncol = res)
maxima <- count_local_maxima(z_mat, x_coords, y_coords)

dx <- diff(x_coords)[1]
dy <- diff(y_coords)[1]

# Gradient components
dz_dx <- (z_mat[, c(2:res, res)] - z_mat[, c(1, 1:(res - 1))]) / (2 * dx)
dz_dy <- (z_mat[c(2:res, res), ] - z_mat[c(1, 1:(res - 1)), ]) / (2 * dy)


coord_to_index <- function(x, y) {
  ix <- which.min(abs(x_coords - x))
  iy <- which.min(abs(y_coords - y))
  c(ix, iy)
}

# Calculate path of steepest ascent
steepest_ascent <- function(x0,
                            y0,
                            step_size = 0.01,
                            n_steps = 2000) {
  path <- matrix(NA, n_steps, 2)
  path[1, ] <- c(x0, y0)
  
  for (i in 2:n_steps) {
    idx <- coord_to_index(path[i - 1, 1], path[i - 1, 2])
    ix <- idx[1]
    iy <- idx[2]
    
    gx <- dz_dx[iy, ix]
    gy <- dz_dy[iy, ix]
    
    grad_norm <- sqrt(gx^2 + gy^2)
    if (grad_norm == 0)
      break
    
    # Move uphill
    path[i, 1] <- path[i - 1, 1] + step_size * gx / grad_norm
    path[i, 2] <- path[i - 1, 2] + step_size * gy / grad_norm
  }
  
  path_df <- as.data.frame(path)
  colnames(path_df) <- c("x", "y")
  path_df
}


ascent_path <- steepest_ascent(6.4, 1.12)


p_A <- ggplot(noise_df, aes(x = x, y = y, fill = z)) +
  geom_raster() +
  
  geom_point(
    data = filter(results, motif != 0, rep == 1, motif == 1 |
                    motif == 5),
    aes(x = mean_x1, y = mean_x2, colour = time),
    #color = "red",
    size = 1,
    inherit.aes = FALSE
  ) +
  geom_path(
    data = ascent_path,
    aes(x = y, y = x),
    color = "red",
    linewidth = 1,
    inherit.aes = FALSE,
  ) +
  
  #coord_equal() +
  scale_fill_viridis_c() +
  facet_wrap(~ motif, labeller = labeller(
    motif = function(x)
      paste("Motif", x)
  )) +
  theme_minimal() +
  annotate(
    "point",
    x = maxima$x[1:2],
    y = maxima$y[1:2],
    shape = 4,
    size = 1,
    color = "red"
  ) +
  labs(
    title = "",
    fill = "Fitness",
    color = "Time",
    x = expression("Phenotype" ~ x[1]),
    y = expression("Phenotype" ~ x[2])
    
    
  )
p_A
#ggsave("motif_div-trajectories.png", width = 7, height = 5, bg = "white")

# Plot B
cov_ellipse <- function(mu, Sigma, level = 0.95, n = 100) {
  r <- sqrt(qchisq(level, df = 2))
  theta <- seq(0, 2 * pi, length.out = n)
  circle <- cbind(cos(theta), sin(theta))
  
  ell <- r * circle %*% chol(Sigma)
  tibble(x = ell[, 1] + mu[1], y = ell[, 2] + mu[2])
}


ellipse_df <- results2 %>%
  filter(motif > 0, time %% interval == 0) %>%
  rowwise() %>%
  mutate(ellipse = list({
    G_matrix <- matrix(c(Gvar1_, Gcov_, Gcov_, Gvar2_),
                    nrow = 2,
                    byrow = TRUE)
    
    if (any(eigen(G_matrix, symmetric = TRUE)$values <= 0)) {
      return(NULL)
    }
    
    cov_ellipse(mu = c(mean_x1, mean_x2), Sigma = 2*G_matrix)
  })) %>%
  unnest(ellipse)



Sigma <- mut_rate/2 * matrix(c(2 * mut_dist, 0, 0, 2 * mut_dist), nrow = 2)


data_for_M <- results %>%
  filter(motif > 0, time %% interval == 0 | time == 1) %>%
  group_by(motif, rep, time) %>%
  summarise(x1 = first(mean_x1),
            x2 = first(mean_x2),
            .groups = "drop")


calculate_M_ellipses_df <- function(data, Sigma, level = 0.95) {
  data %>%
    mutate(ellipse = pmap(list(motif, x1, x2), function(j, x, y) {
      J_fun <- J_matrix[[j]]
      M <- J_fun(x, y, th = 5) %*% (Sigma) %*% t(J_fun(x, y, th = 5))
      
      if (any(is.na(M)) ||
          any(eigen(M, symmetric = TRUE)$values <= 0)) {
        return(NULL)
      }
      
      mat_to_ellipse(M, centre = c(x, y), level)
    })) %>%
    unnest(ellipse, keep_empty = TRUE)
}

M_matrices <- calculate_M_ellipses_df(data_for_M, 2*Sigma)

p_B <- ggplot(filter(results, motif == 1 |
                       motif == 5, rep == 1),
              aes(mean_x1, mean_x2, group = rep)) +
  geom_path() +
  geom_path(
    data = filter(M_matrices, rep == 1, motif == 1 | motif == 5),
    aes(ex, ey, group = interaction(rep, time)),
    inherit.aes = FALSE,
    colour = "#00BFC4",
    alpha = 0.5
  ) +
  geom_path(
    data = filter(ellipse_df, rep == 1, motif == 1 | motif == 5, ),
    aes(x, y, group = interaction(rep, time)),
    inherit.aes = FALSE,
    colour = "#F8766D",
    alpha = 0.5
  ) +
  
  facet_wrap(~ motif) +
  coord_equal() +
  theme_minimal() + xlab(expression(x[1])) + ylab(expression(x[2])) + theme(strip.text = element_blank())

p_B
############

matrix_data <- filter(results2, motif == 1 | motif == 5, rep == 1)
matrix_data$Gvar1 <- matrix_data$Gvar1_
matrix_data$Gvar2 <- matrix_data$Gvar2_
matrix_data$Gcov <- matrix_data$Gcov_

matrix_data <- select(matrix_data, -Gvar1_, -Gvar2_, -Gcov_)

matrix_data <- matrix_data %>%
  rowwise() %>%
  mutate(
    Mvar1 = (J_matrix[[motif]](mean_x1, mean_x2, th = 5) %*% (Sigma) %*% t(J_matrix[[motif]](
      mean_x1, mean_x2, th = 5
    )))[1, 1],
    Mvar2 = (J_matrix[[motif]](mean_x1, mean_x2, th = 5) %*% (Sigma) %*% t(J_matrix[[motif]](
      mean_x1, mean_x2, th = 5
    )))[2, 2],
    Mcov = (J_matrix[[motif]](mean_x1, mean_x2, th = 5) %*% (Sigma) %*% t(J_matrix[[motif]](
      mean_x1, mean_x2, th = 5
    )))[1, 2],
    
  )


matrix_data_renamed <- matrix_data %>%
  pivot_longer(
    cols = c(Mvar1, Mvar2, Mcov, Gvar1, Gvar2, Gcov),
    names_to = c("matrix", ".value"),
    names_pattern = "([MG])(.*)"
  )

matrix_data_renamed <- matrix_data_renamed %>%
  mutate(
    disc = sqrt((var1 - var2)^2 + 4 * cov^2),
    lambda_max = (var1 + var2 + disc) / 2,
    lambda_min = (var1 + var2 - disc) / 2,
    ecc = sqrt(1 - lambda_min / lambda_max)
  )

my_scale <- scale_color_discrete(name = NULL, labels = c("G-matrix", "M-matrix"))

p_C1 <- ggplot(filter(matrix_data_renamed, matrix == "G"), aes(time, ecc, colour = matrix)) +
  geom_line() + my_scale +
  facet_wrap(~ motif, labeller = labeller(
    motif = function(x)
      paste("Motif", x)
  )) + ylab("Eccentricity") + theme_minimal() + theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )


matrix_results <- matrix_data  %>%
  rowwise() %>%
  mutate(
    # Eccentricity
    M_ecc = sqrt(1 - ((
      Mvar1 + Mvar2 - sqrt((Mvar1 - Mvar2)^2 + 4 * Mcov^2)
    ) / 2) /
      ((
        Mvar1 + Mvar2 + sqrt((Mvar1 - Mvar2)^2 + 4 * Mcov^2)
      ) / 2)),
    
    G_ecc = sqrt(1 - ((
      Gvar1 + Gvar2 - sqrt((Gvar1 - Gvar2)^2 + 4 * Gcov^2)
    ) / 2) /
      ((
        Gvar1 + Gvar2 + sqrt((Gvar1 - Gvar2)^2 + 4 * Gcov^2)
      ) / 2)),
    
    # Major axis angle (radians)
    M_angle = ((0.5 * atan2(2 * Mcov, Mvar1 - Mvar2) * 180 / pi) + 180) %% 180,
    G_angle = (0.5 * atan2(2 * Gcov, Gvar1 - Gvar2) * 180 / pi) + 180 %% 180,
    
    # Angle difference (degrees)
    angle_diff = abs(M_angle - G_angle) * 180 / pi
  ) %>%
  ungroup()


p_C2 <- ggplot(matrix_results) +
  geom_line(aes(time, G_angle, colour = "G-matrix")) +
  geom_line(aes(time, M_angle, colour = "M-matrix")) +
  facet_wrap(~ motif) + labs(y = "Angle", colour = "") +  scale_color_discrete(guide = "none") +
  scale_y_continuous(breaks = c(-50, 0, 50),
                     labels = c("-50°", "0°", "50°")) +
  theme_minimal() + theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + theme(strip.text = element_blank())




matrix_results <- matrix_results %>%
  mutate(
    # generalized variance = determinant
    G_gen_var = Gvar1 * Gvar2 - Gcov^2,
    M_gen_var = (Mvar1 * Mvar2 - Mcov^2)
  )

matrix_results <- matrix_results %>%
  mutate(G_total_var = Gvar1 + Gvar2,
         M_total_var = 0.001 * (Mvar1 + Mvar2))

p_C3 <- ggplot(matrix_results) +
  geom_line(aes(time, G_gen_var, colour = "G-matrix")) +
  geom_line(aes(time, M_gen_var, colour = "M-matrix")) +
  facet_wrap(~ motif) + labs(x = "Time", y = "Determinant", colour = "") +  scale_color_discrete(guide = "none") +
  theme_minimal() + theme(strip.text = element_blank()) + xlim(c(0, 10000))




library(patchwork)

p_A <- p_A + labs(tag = "A")
p_B <- p_B + labs(tag = "B")
p_C1 <- p_C1 + labs(tag = "C")


left_column <- p_A / p_B + plot_layout()
right_column <- (p_C1 / p_C2 / p_C3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


combined_patchwork <- (combined_plot | right_column) +
  plot_layout(widths = c(1.5, 1))


combined_patchwork

ggsave("fig6_new.png", width = 8.3, height = 4.65, dpi=300)
ggsave("fig6_new.svg", width = 8.3, height = 4.65, dpi=300)

####################################

library(grid)
library(gtable)
library(cowplot)
library(patchwork)

xlims <- c(0, 12)
ylims <- c(0, 8)

common_theme <- theme_minimal() +
  theme(
    strip.text = element_text(),
    strip.background = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(3, "mm"),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

# --- 4. Rebuild plots with identical structure --------------------------------

p_A2 <- ggplot(noise_df, aes(x = x, y = y, fill = z)) +
  geom_raster() +
  geom_point(
    data = filter(results, motif != 0, rep == 1, motif %in% c(1, 5)),
    aes(x = mean_x1, y = mean_x2, colour = time),
    size = 0.5,
    inherit.aes = FALSE
  ) +
  # geom_path(
  #   data = ascent_path,
  #   aes(x = y, y = x),
  #   color = "red",
  #   linewidth = 1,
  #   inherit.aes = FALSE
  # ) +
  annotate(
    "point",
    x = maxima$x[1:2],
    y = maxima$y[1:2],
    shape = 4,
    size = 0.75,
    stroke = 0.75,
    colour = "red"
  ) +
  scale_fill_viridis_c() +
  scale_x_continuous(limits = xlims, expand = c(0, 0)) +
  scale_y_continuous(limits = ylims, expand = c(0, 0)) +
  facet_wrap(~ motif,
             labeller = labeller(
               motif = function(x)
                 paste("Motif", x)
             ),
             scales = "fixed") +
  coord_equal() +
  labs(x = expression("Phenotype" ~ x[1]), y = expression("Phenotype" ~ x[2])) +
  common_theme + labs(tag = "A")



p_B2 <- ggplot(filter(results, motif %in% c(1, 5), rep == 1),
               aes(mean_x1, mean_x2, group = rep)) +
  geom_path(linewidth = 0.3) +
  #geom_path(
  #  data = filter(M_matrices, rep == 1, motif %in% c(1, 5)),
  #  aes(ex, ey, group = interaction(rep, time)),
  #  colour = "#00BFC4",
  #  alpha = 0.5,
  #  inherit.aes = FALSE
  #) +
  geom_path(
    data = filter(ellipse_df, rep == 1, motif %in% c(1, 5)),
    aes(x, y, group = interaction(rep, time)),
    colour = "#F8766D",
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  annotate(
    "point",
    x = maxima$x[1:2],
    y = maxima$y[1:2],
    shape = 4,
    size = 0.75,
    stroke = 0.75,
    colour = "red"
  ) +
  scale_x_continuous(limits = xlims, expand = c(0, 0)) +
  scale_y_continuous(limits = ylims, expand = c(0, 0)) +
  facet_wrap(~ motif,
             labeller = labeller(
               motif = function(x)
                 ""
             ),
             scales = "fixed") +
  coord_equal() +
  labs(x = expression("Phenotype" ~ x[1]), y = expression("Phenotype" ~ x[2])) +
  common_theme + labs(tag = "B")

# --- 5. Convert to grobs and FORCE alignment ----------------------------------

g1 <- ggplotGrob(p_A2)
g2 <- ggplotGrob(p_B2)

# Force identical widths/heights (this is the key step)
g1$widths  <- unit.pmax(g1$widths, g2$widths)
g2$widths  <- unit.pmax(g1$widths, g2$widths)

g1$heights <- unit.pmax(g1$heights, g2$heights)
g2$heights <- unit.pmax(g1$heights, g2$heights)

# # --- 6. Draw final aligned plot -----------------------------------------------
# label_A <- textGrob("A", x = unit(-0.1, "npc"), y = unit(1.1, "npc"),
#                     just = c("left", "top"),
#                     gp = gpar(fontsize = 16, fontface = "bold"))
# label_A <- grobTree(
#   rectGrob(gp = gpar(fill = NA, col = NA)),
#   textGrob("A", gp = gpar(fontsize = 16, fontface = "bold"))
# )
# label_B <- textGrob("B", x = unit(0, "npc"), y = unit(1, "npc"),
#                     just = c("left", "top"),
#                     gp = gpar(fontsize = 16, fontface = "bold"))
#
# label_A <- grobTree(
#   textGrob("A", gp = gpar(fontsize = 16, fontface = "bold")),
#   cl = "off"
# )
#
# label_B <- grobTree(
#   textGrob("B", gp = gpar(fontsize = 16, fontface = "bold")),
#   cl = "off"
# )
#
# g1 <- gtable_add_grob(g1, label_A,
#                       t = 1, l = 1,
#                       z = Inf)
#
# g2 <- gtable_add_grob(g2, label_B,
#                       t = 1, l = 1,
#                       z = Inf)
#
# g1 <- gtable_add_grob(g1, label_A, t = 1, l = 1, z = Inf)
# g2 <- gtable_add_grob(g2, label_B, t = 1, l = 1, z = Inf)
#


grid.newpage()
#combined <- grid.draw(rbind(g1, g2))
combined <- rbind(g1, g2)
combined_plot <- wrap_elements(full = combined)





################################

# G-matrix determination, hether and hohenlohe


getG <- function(G.object,
                 useAll = "TRUE",
                 totalMales = Pop1.sires,
                 totalFemales = Pop1.dams,
                 Bayes = TRUE) {
  X = G.object$IELs
  motif = G.object$motif
  
  if (G.object$n.sires > design.limit.sires) {
    n.sires <- design.limit.sires
  } else {
    n.sires <- G.object$n.sires
  }
  if (G.object$n.dams > design.limit.dams) {
    n.dams <- design.limit.dams
  } else {
    n.dams <- G.object$n.dams
  }
  parents <- G.object$IELs
  
  
  P1 <- rep(sample(1:totalMales, n.sires, replace = FALSE), each = n.dams)
  P2 <- sample((totalMales + 1):(totalMales + totalFemales),
               n.dams * n.sires,
               replace = FALSE)
  
  Ped <- cbind(1:length(P1) + 5000000, P1, P2)
  
  # SIMULATE MATING
  X.offspring <- t(sapply(1:length(P1), function(x, ...) {
    to.mate <- Ped[x, 2:3]
    
    gametes1 <- rbind(as.numeric(X[to.mate[[1]], ][c(1, 3)]), as.numeric(X[to.mate[[1]], ][c(2, 4)]))
    g1 <- apply(gametes1, 2, function(x) {
      tmp <- sample(x, size = 1)
    })
    
    
    gametes2 <- rbind(as.numeric(X[to.mate[[2]], ][c(1, 3)]), as.numeric(X[to.mate[[2]], ][c(2, 4)]))
    g2 <- apply(gametes2, 2, function(x) {
      tmp <- sample(x, size = 1)
    })
    
    return(c(g1[1], g2[1], g1[2], g2[2]))
    
  }))
  
  Ped <- ped <- rbind(cbind(c(unique(P1), unique(P2)), rep(NA, n.sires +
                                                             (n.sires * n.dams)), rep(NA, n.sires + (n.sires * n.dams))), Ped)
  
  colnames(Ped) = c("animal", "FATHER", "DAM")
  
  Ped <- inverseA(Ped)
  Ped <- Ped$Ainv
  head(Ped[1:5, 1:5])
  
  n.traits = 2
  
  RENTS = GenoPhenoMapping(parents[c(unique(P1), unique(P2)), ],
                           motif = motif,
                           plot = FALSE,
                           theta = theta)
  OFFSPRING = GenoPhenoMapping(X.offspring,
                               motif = motif,
                               plot = FALSE,
                               theta = theta)
  
  Data <- data.frame(ped, rbind(RENTS, OFFSPRING))
  colnames(Data) = c("animal", "sire", "dam", "p1", "p2")
  
  
  # Choice of Priors # # # #
  phen.var = cov(Data[, -c(1:3)])
  
  prior2.8 <- list(G = list(G1 = list(
    V = phen.var,
    n = 2,
    alpha.mu = c(0, 0),
    alpha.V = diag(2) * 1000
  )), R = list(V = diag(2) * 1000, n = (0.002)))
  
  
  # # # # # # # # # # # # # #
  
  # Call MCMCglmm
  model <- MCMCglmm(
    cbind(p1, p2) ~ trait - 1,
    random =  ~ us(trait):animal,
    ,
    rcov =  ~ us(trait):units,
    family = c("gaussian", "gaussian"),
    ginverse = list(animal = Ped),
    data = Data,
    burnin = 1000,
    thin = 25,
    nitt = 13000,
    prior = prior2.8,
    verbose = T,
    pr = TRUE
  ) # set pr==TRUE if you want BVs
  
  
  # Summarize
  if (useAll == "TRUE") {
    G = matrix(
      c(
        posterior.mode(model$VCV[, "p1:p1.animal"]),
        posterior.mode(model$VCV[, "p1:p2.animal"]),
        posterior.mode(model$VCV[, "p1:p2.animal"]),
        posterior.mode(model$VCV[, "p2:p2.animal"])
      ),
      2,
      2,
      byrow = TRUE
    )
    
    P = matrix(
      c(
        posterior.mode(model$VCV[, "p1:p1.animal"] + model$VCV[, "p1:p1.units"]),
        posterior.mode(model$VCV[, "p1:p2.animal"] + model$VCV[, "p1:p2.units"]),
        posterior.mode(model$VCV[, "p1:p2.animal"] + model$VCV[, "p1:p2.units"]),
        posterior.mode(model$VCV[, "p2:p2.animal"] + model$VCV[, "p2:p2.units"])
      ),
      2,
      2,
      byrow = TRUE
    )
    
    R = matrix(
      c(
        posterior.mode(model$VCV[, "p1:p1.units"]),
        posterior.mode(model$VCV[, "p1:p2.units"]),
        posterior.mode(model$VCV[, "p1:p2.units"]),
        posterior.mode(model$VCV[, "p2:p2.units"])
      ),
      2,
      2,
      byrow = TRUE
    )
    
    BV1 <- sapply(ped[, 1], function(x) {
      return(posterior.mode(model$Sol[, paste("animal.p1.animal.", x, sep = "")]))
    })
    BV2 <- sapply(ped[, 1], function(x) {
      return(posterior.mode(model$Sol[, paste("animal.p2.animal.", x, sep = "")]))
    })
    
    h2_p1 = posterior.mode(model$VCV[, "p1:p1.animal"] / (model$VCV[, "p1:p1.animal"] + model$VCV[, "p1:p1.units"]))
    h2_p2 = posterior.mode(model$VCV[, "p2:p2.animal"] / (model$VCV[, "p2:p2.animal"] + model$VCV[, "p2:p2.units"]))
    
    
    
    return(list(
      G = G,
      P = P,
      R = R,
      model = model,
      BVs = cbind(BV1, BV2),
      h2_p1 = h2_p1,
      h2_p2 = h2_p2
    ))
  }
  
}
