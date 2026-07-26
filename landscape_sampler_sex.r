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
G <- 20000 # to check convergence for low freq
mut_rate <- 0.05
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
  
  genotype <- rep(motifsPtoG[[m]](x1, x2), each=2)/2
  
  pop <- matrix(rep(genotype, pop_size), nrow = pop_size, byrow = TRUE)
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
      pop[cbind(parents[,1], sample(1:2, nrow(pop), TRUE))],
      pop[cbind(parents[,2], sample(1:2, nrow(pop), TRUE))],
      pop[cbind(parents[,1], sample(3:4, nrow(pop), TRUE))],
      pop[cbind(parents[,2], sample(3:4, nrow(pop), TRUE))]
    )
    
    # Mutation
    mutations <- rbern(pop_size * 4, mut_rate/2)
    mut_indices <- which(as.logical(mutations)) - 1
    
    ind <- floor(mut_indices / 4) + 1
    loc <- (mut_indices %% 4) + 1
    
    temp_pop <- pop
    
    pop[cbind(ind, loc)] <- pmax(0.0, pop[cbind(ind, loc)] + rnorm(length(ind), mean = 0, sd = mut_dist))
    
    phenotypes <- Map(motifsGtoP[[m]], pop[,1] + pop[,2], pop[,3] + pop[,4])
    
    if(m==5){
      neg_inds <- which(sapply(phenotypes, function(x) any(x < 0)))
      pop[neg_inds, ] <- temp_pop[neg_inds, ]
      phenotypes <- Map(motifsGtoP[[m]], pop[,1] + pop[,2], pop[,3] + pop[,4])
    }
    
    P <- do.call(rbind, phenotypes)
    
    trait1 <- P[,1]
    trait2 <- P[,2]
    
    fitness <- ((gen_perlin(trait1, trait2, frequency = my_freq, seed = my_seed) + 1) / 2)
    
    
    
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
my_frequencies <- c(0.3, 0.6, 1.8)
my_frequencies <- c(1.2)



plan(multisession)  
options(future.rng.onMisuse = "error")
t0<-Sys.time()

for (freq in my_frequencies) {
  for (lan in 1:20) {
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
      xmin = maxima$x - zone,
      xmax = maxima$x + zone,
      ymin = maxima$y - zone,
      ymax = maxima$y + zone
    )
    
    # grid screen
    targets1 <- rep(seq(2, 10, length.out= 9), 9) #41
    targets2 <- rep(seq(2, 10, length.out= 9), each = 9) #41
    
    
    for (m in 1:6) {
      print(Sys.time() - t0)
      t0<-Sys.time()
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
final_results_12 <- final_results
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

wide2 <- wide %>% filter(motif1 != 0, motif6 != 0, motif2!=0, motif3!=0, motif4!=0, motif5!=0)


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
    motif1 = "Motif 1",
    motif2 = "Motif 2",
    motif3 = "Motif 3",
    motif4 = "Motif 4",
    motif5 = "Motif 5",
    motif6 = "Motif 6"
  )))+ xlab("Landscape complexity") + ylab("Outcome divergence")+
  theme_minimal()+theme(legend.position = "none")

ggsave("outcomes2.png", width = 7, height = 5)


ggplot(filter(pairwise_by_landscape, motif_2 == "motif6" ), aes(x= frequency, y = 1 - prop_same, colour = motif_1) ) + geom_point() + 
  facet_wrap(~motif_1, labeller = as_labeller(c(
    motif1 = "Motif 1",
    motif2 = "Motif 2",
    motif3 = "Motif 3",
    motif4 = "Motif 4",
    motif5 = "Motif 5",
    motif6 = "Motif 6"
  )))+ xlab("Landscape complexity") + ylab("Outcome divergence")+
  theme_minimal()+theme(legend.position = "none")

p1 <- ggplot(filter(pairwise_by_landscape, motif_2 == "motif6" ), aes(x = as.factor(frequency),
                                                                      y = 1 - prop_same,
                                                                      color = motif_1)) +
  geom_point(position = position_dodge(width = 0.5)) +scale_color_discrete(
    name = "",
    labels = c("Motif 1", "Motif 2", "Motif 3", "Motif 4", "Motif 5")
  )+
  labs(x = "Landscape complexity",
       y = "Outcome divergence with respect to motif 6")+theme_bw()

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
    offset = 0.2,                # how far to shift outward
    x_label = x + offset * dx,    # node position + outward shift
    y_label = y + offset * dy
  )

layout$y_label[4] <- layout$y_label[4] + 0.3
layout$y_label[4] <- layout$y_label[4] - 0.1
layout$x_label[4] <- layout$x_label[4] - 0.2
layout$x_label[5] <- layout$x_label[5] - 0.2
layout$y_label[5] <- layout$y_label[5] + 0.1
layout$x_label[1] <- layout$x_label[1] - 0.1
layout$y_label[1] <- layout$y_label[1] + 0.1
layout$y_label[5] <- layout$y_label[5] - 0.1
layout$y_label[2] <- layout$y_label[2] - 0.15
layout$y_label[3] <- layout$y_label[3] + 0.05
layout$x_label[6] <- layout$x_label[6] - 0.05
layout$y_label[6] <- layout$y_label[6] - 0.05
layout$y[3] <- layout$y[3] - 0.1
layout$x[6] <- - layout$x[6]
layout$x[6] <- - layout$x[6]
layout$x[6] <- - layout$x[6]

p2 <- ggraph(layout) +
  geom_edge_link(aes(width = prop_same, alpha = prop_same)) +
  geom_node_point(size = 5) +
  geom_text(aes(x = x_label, y = y_label, label = recode(name,
                                                         "motif1" = "Motif 1",
                                                         "motif2" = "Motif 2",
                                                         "motif3" = "Motif 3",
                                                         "motif4" = "Motif 4",
                                                         "motif5" = "Motif 5",
                                                         "motif6" = "Motif 6")),
            size = 4) +
  scale_edge_width(range = c(0.5, 3), name = "Same outcomes") +
  scale_edge_alpha(range = c(0.3, 1),
                   name = "Same outcomes") +
  theme_void() +coord_equal()+
  labs(title = "")+coord_cartesian(clip = "off")+
  theme(
    plot.margin = margin(20, 40, 20, 40)
  )

p2

ggarrange(p1, p2, labels = c("A", "B"), widths = c(1.2,1))


ggsave("fig5_2.png", width = 11, height = 5, bg = "white", dpi = 600)



#################
# variable mutation rates
# Run replicate function, for parallel-threading
run_rep <- function(m, x1, x2, my_freq, my_seed, maxima_, mut_dist1) {
  
  genotype <- rep(motifsPtoG[[m]](x1, x2), each=2)/2
  
  pop <- matrix(rep(genotype, pop_size), nrow = pop_size, byrow = TRUE)
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
      pop[cbind(parents[,1], sample(1:2, nrow(pop), TRUE))],
      pop[cbind(parents[,2], sample(1:2, nrow(pop), TRUE))],
      pop[cbind(parents[,1], sample(3:4, nrow(pop), TRUE))],
      pop[cbind(parents[,2], sample(3:4, nrow(pop), TRUE))]
    )
    
    # Mutation
    mutations <- rbern(pop_size * 4, mut_rate/2)
    mut_indices <- which(as.logical(mutations)) - 1
    
    ind <- floor(mut_indices / 4) + 1
    loc <- (mut_indices %% 4) + 1
    
    temp_pop <- pop
    
    pop[cbind(ind, loc)] <- pmax(0.0, pop[cbind(ind, loc)] + rnorm(length(ind), mean = 0, sd = mut_dist1))
    
    phenotypes <- Map(motifsGtoP[[m]], pop[,1] + pop[,2], pop[,3] + pop[,4])
    
    if(m==5){
      neg_inds <- which(sapply(phenotypes, function(x) any(x < 0)))
      pop[neg_inds, ] <- temp_pop[neg_inds, ]
      phenotypes <- Map(motifsGtoP[[m]], pop[,1] + pop[,2], pop[,3] + pop[,4])
    }
    
    P <- do.call(rbind, phenotypes)
    
    trait1 <- P[,1]
    trait2 <- P[,2]
    
    fitness <- ((gen_perlin(trait1, trait2, frequency = my_freq, seed = my_seed) + 1) / 2)
    
    
    
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
