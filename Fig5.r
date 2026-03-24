#######
# Code for fig. 5: Run simulations for a number of landscapes of a given frequency,
# compare the outcomes
#######


### Libraries
library(tidyverse)
library(deSolve)
library(Rlab)
library(gtools)
library(ambient)
library(future)
library(future.apply)
library(rootSolve)

#Parameters
pop_size = 1000
G <- 10000
mut_rate <- 0.01
mut_dist <- 0.2 # not default parameter, rerun?
threshold <- 500


count_local_maxima <- function(z, x, y) {
  nrow_z <- nrow(z)
  ncol_z <- ncol(z)
  
  z_pad <- matrix(-Inf, nrow = nrow_z + 2, ncol = ncol_z + 2)
  z_pad[2:(nrow_z + 1), 2:(ncol_z + 1)] <- z
  
  maxima <- matrix(FALSE, nrow = nrow_z, ncol = ncol_z)
  
  for (i in 2:(nrow_z + 1)) {
    for (j in 2:(ncol_z + 1)) {
      neighborhood <- z_pad[(i-1):(i+1), (j-1):(j+1)]
      maxima[i-1, j-1] <- z_pad[i, j] == max(neighborhood)
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

# Run replicate function, for parallel-threading
run_rep <- function(m, x1, x2, my_freq, my_seed, maxima_) {

  genotype <- motifsPtoG[[m]](x1, x2)
  
  pop <- replicate(pop_size, genotype, simplify = FALSE)
  phenotypes <- replicate(pop_size, c(x1, x2), simplify = FALSE)
  fitness <-   rep((gen_perlin(
    x = x1 + 0.5,
    y = x2 + 0.5,
    frequency = my_freq,
    seed = my_seed
  ) + 1) / 2, pop_size)
  
  
  results <- data.frame(
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
      
      if(m == 5){
        temp_ind <- pop[[ind]]
        temp_ind[loc]  <- max(0.0, pop[[ind]][loc] + rnorm(1, 0, mut_dist))
        temp_phenotype <- motifsGtoP[[m]](pop[[ind]][[1]], pop[[ind]][[2]])
        if(all(temp_phenotype > 0)){
          pop[[ind]] <- temp_ind
          phenotypes[[ind]] <- temp_phenotype
        }
      }
      else{
        pop[[ind]][loc] <- max(0.0, pop[[ind]][loc] + rnorm(1, 0, mut_dist))
        
        phenotypes[[ind]] <- motifsGtoP[[m]](pop[[ind]][[1]], pop[[ind]][[2]])
      }


      #Updating fitness
      fitness[[ind]] <- (
        gen_perlin(
          x = phenotypes[[ind]][1] + 0.5,
          y = phenotypes[[ind]][2] + 0.5,
          frequency = my_freq,
          seed = my_seed
        ) + 1
      ) / 2
    }
    
    
    #Conditional stopping
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
        results <- bind_rows(results,
                             tibble(
                               x1 = x1, 
                               x2 = x2,
                               peak_reached = first(which(condition)),
                               time    = g,
                               mean_X1 = mean(df$x1),
                               mean_X2 = mean(df$x2)
                               
                             ))
        
        break
      }
    }
    # Final stop
    if (g == G) {
      df <- as.data.frame(do.call(rbind, phenotypes))
      names(df) <- c("x1", "x2")
      
      results <- bind_rows(results,
                           tibble(
                             x1 = x1,
                             x2 = x2,
                             peak_reached = 0,
                             time    = g,
                             mean_X1 = mean(df$x1),
                             mean_X2 = mean(df$x2)
                             
                           ))
    }
  }
  results
}


# Simulation loop
all_results <- list()
k <- 1
my_frequencies <- c(0.1, 0.3, 0.6, 1.2, 1.8)



plan(multisession)  
options(future.rng.onMisuse = "error")

for (freq in my_frequencies) {
  for (lan in 1:10) {
    print(lan)
    my_seed <- sample(1:10000, 1) # random landscape / seed
    
    res <- 1000 # Resolution (100x100)
    x_coords <- seq(0.5, 30.5, length.out = res) # shifted grid to avoid lattice effects from perlin noise
    y_coords <- seq(0.5, 30.5, length.out = res)
    noise_matrix <- outer(x_coords, y_coords, function(x, y) {
      gen_perlin(x, y, frequency = freq, seed = my_seed)
    })

    # determine local maxima
    maxima <- count_local_maxima(noise_matrix, x_coords, y_coords)

    #zone within which pop has to cross threshold
    zone <- 0.1
    maxima_ <- data.frame(
      xmin = maxima$x - zone - 0.5,
      xmax = maxima$x + zone - 0.5,
      ymin = maxima$y - zone - 0.5,
      ymax = maxima$y + zone - 0.5
    )
    
    # grid screen
    targets1 <- rep(seq(2, 10, length.out= 41), 41)
    targets2 <- rep(seq(2, 10, length.out= 41), each = 41)
    
    
    for (m in 1:6) {
      df_list <- future_lapply(seq_along(targets1), function(i_)
        run_rep(m, targets1[i_], targets2[i_], freq, my_seed, maxima_), future.seed = TRUE)
      
      results_ <- do.call(rbind, df_list)
      
      results_$motif <- m
      results_$landscape <- lan
      results_$frequency <- freq
      all_results[[k]] <- results_
      k <- k + 1
    }
  }
}
Sys.time()

final_results <- do.call(rbind, all_results)

cleaned <- final_results %>%
  group_by(frequency, landscape, motif, x1, x2) %>%
  slice(1) %>%      # keep only the first row per group
  ungroup()

wide <- cleaned %>%
  distinct(frequency, landscape, motif, x1, x2, peak_reached) %>%
  pivot_wider(
    names_from = motif,
    values_from = peak_reached,
    names_prefix = "motif"
  )

wide2 <- wide %>% filter(motif1 != 0, motif6 != 0)


motif_cols <- grep("^motif", names(wide), value = TRUE)

pairwise_by_landscape <- wide2 %>%
  group_by(frequency, landscape) %>%
  group_modify(~ {
    map_dfr(
      combn(motif_cols, 2, simplify = FALSE),
      function(cols) {
        m1 <- cols[1]
        m2 <- cols[2]
        
        .x %>%
          filter(!is.na(.data[[m1]]), !is.na(.data[[m2]])) %>%
          summarise(
            motif_1 = m1,
            motif_2 = m2,
            n_same  = sum(.data[[m1]] == .data[[m2]]),
            n_total = n(),
            prop_same = n_same / n_total
          )
      }
    )
  }) %>%
  ungroup()


plot_df2 <- pairwise_by_landscape %>% group_by(frequency, landscape)%>%
  summarize(
    mean_prop = mean(prop_same)
  ) %>%
  ungroup()

p1 <- ggplot(plot_df2, aes(x = frequency, y = mean_prop))+geom_point() + xlab("Landscape complexity") + ylab("Outcome consistency")+theme_minimal()

ggsave("fig6_outcomes.png", width = 7, height = 5)


plot_df <- pairwise_by_landscape %>%
  mutate(
    comparison = paste(motif_1, motif_2, sep = " – ")
  )



ggplot(filter(pairwise_by_landscape, motif_2 == "motif6" ), aes(x= frequency, y = 1 - prop_same, colour = motif_1) ) + geom_point() + 
  facet_wrap(~motif_1, labeller = as_labeller(c(
    motif1 = "motif 1",
    motif2 = "motif 2",
    motif3 = "motif 3",
    motif4 = "motif 4",
    motif5 = "motif 5",
    motif6 = "motif 6"
    )))+ xlab("Landscape complexity") + ylab("Outcome divergence")+
  theme_minimal()+theme(legend.position = "none")

ggsave("outcomes2.png", width = 7, height = 5)


### # Network comparison, star

library(dplyr)
library(purrr)

motif_cols <- grep("^motif", names(wide2), value = TRUE)

pairwise_results <- map_dfr(
  combn(motif_cols, 2, simplify = FALSE),
  function(cols) {
    m1 <- cols[1]
    m2 <- cols[2]
    
    wide %>%
      filter(!is.na(.data[[m1]]), !is.na(.data[[m2]])) %>%
      summarise(
        motif_1 = m1,
        motif_2 = m2,
        n_same  = sum(.data[[m1]] == .data[[m2]]),
        n_total = n(),
        prop_same = n_same / n_total
      )
  }
)

library(igraph)
library(ggraph)

graph <- graph_from_data_frame(
  pairwise_results,
  directed = FALSE,
  vertices = unique(c(pairwise_results$motif_1, pairwise_results$motif_2))
)


layout <- create_layout(graph, layout = "fr")
layout <- layout %>%
  mutate(
    r = sqrt(x^2 + y^2),          # distance from center
    dx = x / r,                   # outward direction (unit vector)
    dy = y / r,
    offset = 0.15,                # how far to shift outward
    x_label = x + offset * dx,    # node position + outward shift
    y_label = y + offset * dy
  )

layout$y_label[4] <- layout$y_label[4] - 0.1
layout$y_label[1] <- layout$y_label[1] - 0.1
layout$y_label[5] <- layout$y_label[5] - 0.2
layout$y_label[2] <- layout$y_label[2] - 0.05
layout$y_label[3] <- layout$y_label[3] - 0.05
layout$y_label[6] <- layout$y_label[6] - 0.05
layout$y[3] <- layout$y[3] - 0.1
layout$x[6] <- - layout$x[6]
layout$x[6] <- - layout$x[6]
layout$x[6] <- - layout$x[6]

ggraph(layout) +
  geom_edge_link(aes(width = prop_same, alpha = prop_same)) +
  geom_node_point(size = 5) +
  geom_text(aes(x = x_label, y = y_label, label = recode(name,
                                                         "motif1" = "motif 1",
                                                         "motif2" = "motif 2",
                                                         "motif3" = "motif 3",
                                                         "motif4" = "motif 4",
                                                         "motif5" = "motif 5",
                                                         "motif6" = "motif 6")),
            size = 4) +
  scale_edge_width(range = c(0.5, 3), name = "Same outcomes") +
  scale_edge_alpha(range = c(0.3, 1),
                   name = "Same outcomes") +
  theme_void() +coord_equal()+
  labs(title = "")+coord_cartesian(clip = "off")+
  theme(
    plot.margin = margin(20, 40, 20, 40)
  )

ggsave("fig6_2.png", width = 7, height = 5.5, bg = "white")

