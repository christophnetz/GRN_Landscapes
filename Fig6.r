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
  
  x1 <- a1/g
  x2 <- a2 * (x1/(th+x1))/g
  
  return(c(x1, x2))
}

motif_2 <- function(a1, a2, g = 1, th = 5) {
  
  x1 <- a1/g
  x2 <- (a2*(1-((a1/g)/(th+((a1)/(g))))))/(g)
  
  return(c(x1, x2))
}

motif_3 <- function(a1, a2, g = 1, th = 5) {
  
  x1 <- (2*th*a1*sqrt(g))/((-a1+th*g)*sqrt(g) + sqrt(4*a1*a2*g + ((a1 + th*g)^2)*g))
  x2 <- (-(a1+th*g)*sqrt(g) + sqrt(4*a1*a2*g + ((a1 + th*g)^2)*g))/(2*g*sqrt(g))
  
  return(c(x1, x2))
}

motif_4 <- function(a1, a2, g = 1, th = 5) {
  
  x1 <- (a1-a2-g*th+ sqrt(4*a1*g*th + (a1-a2-g*th)^2))/(2*g)
  x2 <- (-a1+a2-g*th + sqrt(4*a1*g*th + (a1-a2-g*th)^2))/(2*g)
  
  return(c(x1, x2))
}

motif_5 <- function(a1, a2, g = 1, th = 5) {
  
  x1 <- (a1*a2-(g^2)*(th^2))/(g*(a2+g*th))
  x2 <- (a1*a2-(g^2)*(th^2))/(g*(a1+g*th))
  
  return(c(x1, x2))
}

motif_6 <- function(a1, a2, g = 1, th = 5) {
  
  x1 <- a1/g
  x2 <- a2/g
  
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
  
  a1 <- g*x1
  a2 <- (g*x2*(th+x1))/(x1)
  return(c(a1, a2))
}

motifPG_2 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- g*x1
  a2 <- (g*x2)/(1-(x1/(th+x1)))
  return(c(a1, a2))
}

motifPG_3 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- (g*x1)/(1-((x2)/(th+x2)))
  a2 <- (g*x2*(th+x1))/(x1)
  
  return(c(a1, a2))
}

motifPG_4 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- (g*x1*(th+x2))/(th)
  a2 <- (g*x2*(th+x1))/(th)
  return(c(a1, a2))
}

motifPG_5 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- (g*x1*(th+x2))/(x2)
  a2 <- (g*x2*(th+x1))/(x1)
  return(c(a1, a2))
}

motifPG_6 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- x1*g
  a2 <- x2*g
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


### Parameters
replicates <- 1
nr_motives <- 2
length_df <- (G)*replicates*nr_motives
results <- data.frame(
  motif  = numeric(length_df),
  rep     = numeric(length_df),
  time    = numeric(length_df),
  mean_x1 = numeric(length_df),
  mean_x2 = numeric(length_df),
  Gvar1   = numeric(length_df),
  Gvar2   = numeric(length_df),
  Gcov = numeric(length_df)
  )

pop_size = 10000
G <- 10000
mut_rate <- 0.01 
mut_dist <- 0.05 


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
for (m in c(1,5)) { # for (m in 1:nr_motives) {
  print(m)

  for (rep in 1:replicates) {
    print(rep)
    
    genotype <- motifsPtoG[[m]](1.12, 6.4)
    
    pop <- replicate(pop_size, genotype, simplify = FALSE)
    phenotypes <- replicate(pop_size, c(1.12, 6.4), simplify = FALSE)
    fitness <-   rep( fitness_fun(((gen_perlin(
      x = 1.12,
      y = 6.4 + 1,
      frequency = 0.05,
      seed = 76
    ) +1) /2)^2, 
    ((gen_perlin(
      x = 1.12 - 1,
      y = 6.4,
      frequency = 0.05,
      seed = 76
    ) +1) /2)^2), pop_size)
    
    
    # generation loop
    for (g in 1:G) {
      
      #sample offspring in a weighted lottery
      inds <- sample(1:pop_size,
                     size = pop_size,
                     replace = TRUE,
                     prob = (fitness - min(fitness) + 0.001))
      
      pop <- pop[inds]
      phenotypes <- phenotypes[inds]
      fitness <- fitness[inds]
      
      # Mutation
      mutations <- rbern(pop_size * 2, mut_factor * mut_rate)
      mut_indices <- which(as.logical(mutations)) - 1
      
      for (mut in mut_indices) {
        ind <- floor(mut / 2) + 1
        loc <- (mut %% 2) + 1
        
        pop[[ind]][loc] <- max(0.0, pop[[ind]][loc] + rnorm(1, 0, mut_dist))
        
        phenotypes[[ind]] <- motifsGtoP[[m]](pop[[ind]][[1]], pop[[ind]][[2]])
        
        #Updating fitness
        fitness[[ind]] <- fitness_fun(
          ((gen_perlin(
            x = phenotypes[[ind]][1],
            y = phenotypes[[ind]][2] + 1,
            frequency = 0.05,
            seed = 76
          )+1)/2)^2, ((gen_perlin(
            x = phenotypes[[ind]][1] - 1,
            y = phenotypes[[ind]][2],
            frequency = 0.05,
            seed = 76
          )+1)/2)^2
        )
      }
      
      output_data <- as.data.frame(do.call(rbind, phenotypes))
      names(output_data) <- c("x1", "x2")
      
      results[(G)*(rep-1) +(G)*(replicates)*(m-1) + g, ] <- c(
        m,
        rep,
        g,
        mean(output_data$x1),
        mean(output_data$x2),
        var(output_data$x1),
        var(output_data$x2),
        cov(output_data$x1, output_data$x2)
        )
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
  gen_perlin(x, y+1, frequency = 0.05, seed = 76)
})

base <- (base + 1)/2

terrain_steep <- base^2

base2 <- outer(x_coords, y_coords, function(x, y) {
  gen_perlin(x-1, y, frequency = 0.05, seed = 76)
})

base2 <- (base2 + 1)/2

terrain_steep2 <- base2^2

softmax2 <- function(a, b, beta = 5) {
  ea <- exp(beta * a)
  eb <- exp(beta * b)
  ea / (ea + eb)
}

w <- softmax2(terrain_steep, terrain_steep2, beta = 100)

terrain <- w * terrain_steep + (1 - w) * terrain_steep2


noise_df <- expand.grid(
  x = x_coords,
  y = y_coords
)

noise_df$z <- as.vector(terrain)


# Recreate matrix from data frame
z_mat <- matrix(noise_df$z, nrow = res, ncol = res)

dx <- diff(x_coords)[1]
dy <- diff(y_coords)[1]

# Gradient components
dz_dx <- (z_mat[, c(2:res, res)] - z_mat[, c(1, 1:(res-1))]) / (2 * dx)
dz_dy <- (z_mat[c(2:res, res), ] - z_mat[c(1, 1:(res-1)), ]) / (2 * dy)


coord_to_index <- function(x, y) {
  ix <- which.min(abs(x_coords - x))
  iy <- which.min(abs(y_coords - y))
  c(ix, iy)
}

# Calculate path of steepest ascent
steepest_ascent <- function(x0, y0, step_size = 0.01, n_steps = 2000) {
  
  path <- matrix(NA, n_steps, 2)
  path[1, ] <- c(x0, y0)
  
  for (i in 2:n_steps) {
    idx <- coord_to_index(path[i-1, 1], path[i-1, 2])
    ix <- idx[1]
    iy <- idx[2]
    
    gx <- dz_dx[iy, ix]
    gy <- dz_dy[iy, ix]
    
    grad_norm <- sqrt(gx^2 + gy^2)
    if (grad_norm == 0) break
    
    # Move uphill
    path[i, 1] <- path[i-1, 1] + step_size * gx / grad_norm
    path[i, 2] <- path[i-1, 2] + step_size * gy / grad_norm
  }
  
  path_df <- as.data.frame(path)
  colnames(path_df) <- c("x", "y")
  path_df
}


ascent_path <- steepest_ascent(6.4, 1.12)


pl_A <- ggplot(noise_df, aes(x = x, y = y, fill = z)) +
  geom_raster() +
  
  geom_point(
    data = filter(results, motif !=0, rep == 1, motif == 1 | motif == 5),
    aes(x = mean_x1, y = mean_x2, colour = time),
    #color = "red",
    size = 1,
    inherit.aes = FALSE
  )+
  geom_path(
    data = ascent_path,
    aes(x = y, y = x),
    color = "red",
    linewidth = 1,
    inherit.aes = FALSE,
  )+
  
  coord_equal() +
  scale_fill_viridis_c() +
  facet_wrap(~motif, labeller=labeller(motif = function(x) paste("motif", x)))+
  theme_minimal() +
  annotate("point", x = c(9.14, 10.34), y = c(1.56, 2.77), shape = 4, size = 1, stroke = 1) +
  labs(
    title = "",
    fill = "fitness",
    x = expression(x[1]),
    y = expression(x[2])
    
    
  )

#ggsave("motif_div-trajectories.png", width = 7, height = 5, bg = "white")

# Plot B
cov_ellipse <- function(mu, Sigma, level = 0.95, n = 100) {
  r <- sqrt(qchisq(level, df = 2))
  theta <- seq(0, 2 * pi, length.out = n)
  circle <- cbind(cos(theta), sin(theta))
  
  ell <- r * circle %*% chol(Sigma)
  tibble(
    x = ell[,1] + mu[1],
    y = ell[,2] + mu[2]
  )
}


ellipse_df <- results %>%
  filter(motif > 0, time %% 200 == 0) %>%
  rowwise() %>%
  mutate(
    ellipse = list({
      Sigma <- matrix(
        c(Gvar1, Gcov,
          Gcov, Gvar2),
        nrow = 2,
        byrow = TRUE
      )
      
      if (any(eigen(Sigma, symmetric = TRUE)$values <= 0)) {
        return(NULL)
      }
      
      cov_ellipse(
        mu = c(mean_x1, mean_x2),
        Sigma = 5 * Sigma
      )
    })
  ) %>%
  unnest(ellipse)



Sigma <- matrix(c(2 * mut_dist, 0, 0, 2 * mut_dist), nrow = 2)


data_for_M <- results %>%
  filter(motif > 0, time %% 200 == 0 | time == 1) %>%
  group_by(motif, rep, time) %>%
  summarise(
    x1 = first(mean_x1),
    x2 = first(mean_x2),
    .groups = "drop"
  )


calculate_M_ellipses_df <- function(data, Sigma, level = 0.3) {
  
  data %>%
    mutate(
      ellipse = pmap(
        list(motif, x1, x2, mut_factor),
        function(j, x, y, m) {
          J_fun <- J_matrix[[j]]
          M <- J_fun(x, y, th = 5) %*% (m * Sigma) %*% t(J_fun(x, y, th = 5))
          
          if (any(is.na(M)) || any(eigen(M, symmetric = TRUE)$values <= 0)) {
            return(NULL)
          }
          
          mat_to_ellipse(M, centre = c(x, y), level)
        }
      )
    ) %>%
    unnest(ellipse, keep_empty = TRUE)
}

M_matrices <- calculate_M_ellipses_df(data_for_M, Sigma)

pl_B <- ggplot(filter(results, motif == 1 | motif == 5, rep == 1), aes(mean_x1, mean_x2, group = rep)) +
  geom_path() +
  geom_path(
    data =filter(M_matrices, rep ==1, motif == 1 | motif == 5),
    aes(ex, ey, group = interaction(rep, time)),
    inherit.aes = FALSE,
    colour = "#00BFC4",
    alpha = 0.5
  ) +  
  geom_path(
    data = filter(ellipse_df, rep == 1, motif == 1 | motif == 5,),
    aes(x, y, group = interaction(rep, time)),
    inherit.aes = FALSE,
    colour = "#F8766D",
    alpha = 0.5
  ) +
  
  facet_wrap(~ motif) +
  coord_equal() +
  theme_minimal()+xlab(expression(x[1])) + ylab(expression(x[2]))+theme(strip.text = element_blank())

############

matrix_data <- filter(results, motif == 1 | motif == 5, rep == 1)

matrix_data <- matrix_data %>% 
  rowwise() %>%
  mutate(
    
    Mvar1 = (J_matrix[[motif]](mean_x1, mean_x2, th = 5) %*% ( Sigma) %*% t(J_matrix[[motif]](mean_x1, mean_x2, th = 5)))[1,1],
    Mvar2 = (J_matrix[[motif]](mean_x1, mean_x2, th = 5) %*% (Sigma) %*% t(J_matrix[[motif]](mean_x1, mean_x2, th = 5)))[2,2],
    Mcov = (J_matrix[[motif]](mean_x1, mean_x2, th = 5) %*% (Sigma) %*% t(J_matrix[[motif]](mean_x1, mean_x2, th = 5)))[1,2],
    
  )


matrix_data_renamed <- matrix_data %>%
  pivot_longer(
    cols = c(Mvar1, Mvar2, Mcov, Gvar1, Gvar2, Gcov),
    names_to = c("matrix", ".value"),
    names_pattern = "([MG])(.*)"
  )

matrix_data_renamed <- matrix_data_renamed %>%
  mutate(
    disc = sqrt((var1 - var2)^2 + 4*cov^2),
    lambda_max = (var1 + var2 + disc)/2,
    lambda_min = (var1 + var2 - disc)/2,
    ecc = sqrt(1 - lambda_min / lambda_max)
  )

my_scale <- scale_color_discrete(
  name = NULL,
  labels = c("G-matrix", "M-matrix")
)

p_C1 <- ggplot(matrix_data_renamed, aes(time, ecc, colour = matrix))+
  geom_line()+ my_scale+ 
  facet_wrap(~motif, labeller = labeller(motif = function(x) paste("motif", x)))+ylab("eccentricity")+theme(axis.title.x = element_blank(),
                                                                                                            axis.text.x  = element_blank(),
                                                                                                            axis.ticks.x = element_blank())

matrix_results <- matrix_data  %>%
  rowwise() %>%
  mutate(
    # Eccentricity
    M_ecc = sqrt(1 - ( (Mvar1 + Mvar2 - sqrt((Mvar1 - Mvar2)^2 + 4*Mcov^2))/2 ) /
                   ( (Mvar1 + Mvar2 + sqrt((Mvar1 - Mvar2)^2 + 4*Mcov^2))/2 )),
    
    G_ecc = sqrt(1 - ( (Gvar1 + Gvar2 - sqrt((Gvar1 - Gvar2)^2 + 4*Gcov^2))/2 ) /
                   ( (Gvar1 + Gvar2 + sqrt((Gvar1 - Gvar2)^2 + 4*Gcov^2))/2 )),
    
    # Major axis angle (radians)
    M_angle = 0.5 * atan2(2*Mcov, Mvar1 - Mvar2),
    G_angle = 0.5 * atan2(2*Gcov, Gvar1 - Gvar2),
    
    # Angle difference (degrees)
    angle_diff = abs(M_angle - G_angle) * 180 / pi
  ) %>%
  ungroup()


p_C2 <- ggplot(matrix_results)+
  geom_line(aes(time, G_angle, colour="G-matrix"))+
  geom_line(aes(time, M_angle, colour="M-matrix"))+
  facet_wrap(~motif)+labs(y = "angle", colour = "")+  scale_color_discrete(guide = "none")+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())+theme(strip.text = element_blank())



matrix_results <- matrix_results %>%
  mutate(
    # generalized variance = determinant
    G_gen_var = Gvar1 * Gvar2 - Gcov^2,
    M_gen_var = 0.001*( Mvar1 * Mvar2 - Mcov^2)
  )

matrix_results <- matrix_results %>%
  mutate(
    G_total_var = Gvar1 + Gvar2,
    M_total_var = 0.001 * (Mvar1 + Mvar2)
  )

p_C3 <- ggplot(matrix_results)+  
  geom_line(aes(time, G_total_var, colour="G-matrix"))+
  geom_line(aes(time, M_total_var, colour="M-matrix"))+
  facet_wrap(~motif)+labs(y = "total variance", colour = "")+  scale_color_discrete(guide = "none")+
  theme(strip.text = element_blank())




library(patchwork)

p_A <- p_A + labs(tag = "A")
p_B <- p_B + labs(tag = "B")
p_C1 <- p_C1 + labs(tag = "C")


left_column <- p_A / p_B
right_column <- (p_C1 / p_C2 / p_C3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_patchwork <- left_column | right_column

combined_patchwork

ggsave("fig6_new.png", width = 8, height = 5)

