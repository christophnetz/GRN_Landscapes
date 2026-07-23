# Figure 2

library(tidyr)
library(purrr)
library(ggplot2)
library(dplyr)
library(ellipse)


### Jacobians
mut_dist <- 0.05

Sigma <- matrix(c(2 * mut_dist, 0, 0, 2 * mut_dist), nrow = 2)

J1 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- g*x1
  a2 <- (g*x2*(th+x1))/(x1)
  
  return(matrix(c(
    1/g,
    a2 * th / (g * (a1 + th)^2), 
    0, 
    a1/(g*(th + a1))
  ), nrow = 2))
}

J2 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- g*x1
  a2 <- (g*x2)/(1-(x1/(th+x1)))
  
  x1 <- a1/g
  x2 <- (a2*(1-((a1/g)/(th+((a1)/(g))))))/(g)
  
  
  return(matrix(c(
    1/g, 
    - a2 * th/(a1 + g * th)^2, 
    0, 
    (1-((a1/g)/(th+((a1)/(g)))))/(g)
  ), nrow = 2))
}

J3 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- (g*x1)/(1-((x2)/(th+x2)))
  a2 <- (g*x2*(th+x1))/(x1)
  
  return(matrix(c((th*(a1 + 2*a2+g*th + sqrt(4 * a1 * a2 +(a1+g*th)^2)))/(2*(a2 + g * th)*sqrt(4 * a1 * a2 + (a1 + g * th)^2)), 
                  ((a1 + 2 * a2 + g * th)/sqrt(4 * a1 * a2 + (a1 + g * th)^2) -1)/(2*g), 
                  -(4 * a1^2 * th)/(sqrt(4 * a1 * a2 + (a1 + g * th)^2)*(-a1+g * th + sqrt(4 * a1 * a2 + (a1 + g * th)^2))^2), 
                  a1/(g * sqrt(4 * a1 * a2 + (a1 + g * th)^2))), nrow = 2))
}

J4 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- (g*x1*(th+x2))/(th)
  a2 <- (g*x2*(th+x1))/(th)
  
  return(matrix(c(
    (1+(a1 - a2 + g * th)/(sqrt(4 * a1 * g * th + (-a1 + a2 + g * th)^2)))/(2*g), 
    (-1+(a1 - a2 + g * th)/(sqrt(4 * a1 * g * th + (-a1 + a2 + g * th)^2)))/(2*g), 
    (-1+(-a1 + a2 + g * th)/(sqrt(4 * a1 * g * th + (-a1 + a2 + g * th)^2)))/(2*g), 
    (1+(-a1 + a2 + g * th)/(sqrt(4 * a1 * g * th + (-a1 + a2 + g * th)^2)))/(2*g)
  ), nrow = 2))
}

J5 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- (g*x1*(th+x2))/(x2)
  a2 <- (g*x2*(th+x1))/(x1)
  
  return(matrix(c(
    a2/(g*(a2 + th)), 
    (th*(a2 + th))/(g*(a1 + th)^2), 
    (th*(a1 + th))/(g*(a2 + th)^2), 
    a1/(g * (a1 + th))
  ), nrow = 2))
}

J6 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- x1*g
  a2 <- x2*g
  
  return(matrix(c(1/g, 0, 0, 1/g), nrow = 2))
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

ellipse_df <- function(x_vals, y_vals, Sigma, J_matrix, level = 0.95) {
  
  expand_grid(
    J_id = seq_along(J_matrix),
    x    = x_vals,
    y    = y_vals
  ) %>%
    mutate(
      ellipse = pmap(
        list(J_id, x, y),
        function(j, x, y) {
          J_fun <- J_matrix[[j]]
          M <- J_fun(x, y, th = 5) %*% Sigma %*% t(J_fun(x, y, th = 5))
          if (any(is.na(M)) || any(eigen(M, symmetric = TRUE)$values <= 0)) {
            return(NULL)
          }
          mat_to_ellipse(M, centre = c(x, y), level)
        }
      )
    ) %>%
    unnest(ellipse, keep_empty = TRUE)
}


ell_df <- ellipse_df(
  x_vals = seq(2, 10, length.out = 5),
  y_vals = seq(2, 10, length.out = 5),
  Sigma  = Sigma,
  J_matrix = J_matrix
)

ell_df <- ell_df %>%
  rename(motif = J_id)

label_motif <- function(value) {
  paste("Motif", value)  
}

sensitivity <- expand_grid(
  motif = 1:6,
  x1    = seq(2, 10, length.out = 5),
  x2    = seq(2, 10, length.out = 5)
)


sensitivity <- sensitivity %>%
  rowwise() %>%
  mutate(
    J = list(J_matrix[[motif]](x1, x2)),
    x1a1 = x1 +  J[1,1],
    x2a1 = x2 +  J[2,1],
    x1a2 = x1 +  J[1,2],
    x2a2 = x2 +  J[2,2]
  ) %>%
  dplyr::select(-J) %>%
  ungroup()


ggplot(ell_df, aes(ex, ey, group = interaction(motif, x, y))) +
  geom_path() +
  geom_segment(data = sensitivity, aes(
    x=x1,
    y=x2,
    xend=x1a1,
    yend = x2a1
  ), arrow = arrow(length = unit(0.07, "cm"), type = "closed"), colour="blue", linewidth = 0.4,   inherit.aes = FALSE)+
  geom_segment(data = sensitivity, aes(
    x=x1,
    y=x2,
    xend=x1a2,
    yend = x2a2
  ), arrow = arrow(length = unit(0.07, "cm"), type = "closed"), colour="red", linewidth = 0.4,   inherit.aes = FALSE)+
  coord_equal() +
  theme_minimal()+
  facet_wrap(~motif, labeller = labeller(motif = label_motif))+xlab(expression("Phenotype" ~ x[1]))+ylab(expression("Phenotype" ~ x[2]))+
  theme(strip.text = element_text(face = "bold"))+
  theme(panel.spacing.x = unit(3, "lines"))  

ggsave("hohenlohe_M-matrix.png", width = 7, height =4, bg = "white")


p1 <- ggplot(filter(ell_df, motif < 3), aes(ex, ey, group = interaction(motif, x, y))) +
  geom_path() +
  geom_segment(data = filter(sensitivity, motif < 3), aes(
    x=x1,
    y=x2,
    xend=x1a1,
    yend = x2a1
  ), arrow = arrow(length = unit(0.07, "cm"), type = "closed"), colour="blue", linewidth = 0.4,   inherit.aes = FALSE)+
  geom_segment(data = filter(sensitivity, motif < 3), aes(
    x=x1,
    y=x2,
    xend=x1a2,
    yend = x2a2
  ), arrow = arrow(length = unit(0.07, "cm"), type = "closed"), colour="red", linewidth = 0.4,   inherit.aes = FALSE)+
  coord_equal() +
  theme_minimal()+
  facet_wrap(~motif, labeller = labeller(motif = label_motif))+xlab(expression("Phenotype" ~ x[1]))+ylab(expression("Phenotype" ~ x[2]))+
  theme(strip.text = element_text(face = "bold"))+
  theme(panel.spacing.x = unit(3, "lines"))  

my_plot <- readRDS("C:/Users/user/Downloads/my_plot.rds")

ggarrange(p1, my_plot)


ggsave("hohenlohe_M-matrix.png", width = 7, height =4, bg = "white")




ggplot(ell_df, aes(ex, ey, group = interaction(motif, x, y))) +
  geom_path() +
  geom_segment(data = sensitivity, aes(
    x=x1,
    y=x2,
    xend=x1a1,
    yend = x2a1
  ), arrow = arrow(length = unit(0.07, "cm"), type = "closed"), colour="blue", linewidth = 0.4,   inherit.aes = FALSE)+
  geom_segment(data = sensitivity, aes(
    x=x1,
    y=x2,
    xend=x1a2,
    yend = x2a2
  ), arrow = arrow(length = unit(0.07, "cm"), type = "closed"), colour="red", linewidth = 0.4,   inherit.aes = FALSE)+
  coord_equal() +
  theme_minimal()+
  facet_wrap(~motif, labeller = labeller(motif = label_motif))+xlab(expression(x[1]))+ylab(expression(x[2]))+
  theme(strip.text = element_blank(), panel.spacing.y = unit(2, "cm") )+
  theme(panel.spacing.x = unit(3, "lines"))  

ggsave("hohenlohe_M-matrix2.png", width = 7, height =5, bg = "white")




