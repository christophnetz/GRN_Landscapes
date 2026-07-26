library(ggpubr)

head(final_results)

final_results2 <- final_results2 %>%
  mutate(
    divergence = sqrt((mean_X1 - x1)^2 + (mean_X2 - x2)^2),
    div_x1= (mean_X1 - x1),
    div_x2= (mean_X2 - x2)
    
  )


ggplot(final_results2, aes(div_x1, div_x2))+
  geom_point(alpha=0.1, size = 0.3)+
  facet_wrap(~motif)+
  geom_smooth(method = "lm", se = T)+ggtitle("frequency 0.5")+
  stat_regline_equation(label.y = 2.7)+labs(x="divergence x1", y="divergence x2")

#final_results_f1 <- final_results

ggplot(final_results_f1, aes(div_x1, div_x2))+
  geom_point(alpha=0.1, size = 0.3)+
  facet_wrap(~motif)+
  geom_smooth(method = "lm", se = T)+ggtitle("frequency 1.0")+
  stat_regline_equation(label.y = 1.4)+labs(x="divergence x1", y="divergence x2")

###################################
#

# M=
Sigma <- matrix(c(2 * mut_dist, 0, 0, 2 * mut_dist), nrow = 2)

M_fun <- function(motif, x1, x2){
  
  J_matrix[[motif]](x1,x2) %*% Sigma %*% t(J_matrix[[motif]](x1,x2))
  
}



df_f1 <- final_results_f1 %>%
  rowwise() %>%
  mutate(
    e_d = as.numeric(
      t(c(div_x1, div_x2)) %*% M_fun(motif, x1, x2) %*% c(div_x1, div_x2) /
        (t(c(div_x1, div_x2)) %*% c(div_x1, div_x2))
    ),
    freq = 1
  ) %>%
  ungroup() %>% dplyr::select(motif, x1, x2, e_d, freq)

df_f05 <- final_results2 %>%
  rowwise() %>%
  mutate(
    e_d = as.numeric(
      t(c(div_x1, div_x2)) %*% M_fun(motif, x1, x2) %*% c(div_x1, div_x2) /
        (t(c(div_x1, div_x2)) %*% c(div_x1, div_x2))
    ),
    freq = 0.5
  ) %>%
  ungroup() %>% dplyr::select(motif, x1, x2, e_d, freq)

df_freq <- rbind(df_f1, df_f05)

ggplot(df_freq, aes(x=as.factor(freq), y=e_d)) + geom_boxplot()+facet_wrap(~motif)


#######################################
# with landscape sampling
load("~/outcomes-long-run.RData")

head(final_results)

final_results <- final_results %>%
  mutate(
    div_x1= (mean_X1 - x1),
    div_x2= (mean_X2 - x2)
    
  )


ggplot(final_results, aes(div_x1, div_x2))+
  geom_point(alpha=0.1, size = 0.3)+
  facet_grid(frequency~motif)+
  geom_smooth(method = "lm", se = T)+
  stat_regline_equation(label.y = 10)+labs(x="divergence x1", y="divergence x2")


df_f <- final_results %>%
  rowwise() %>%
  mutate(
    e_d = as.numeric(
      t(c(div_x1, div_x2)) %*% M_fun(motif, x1, x2) %*% c(div_x1, div_x2) /
        ((t(c(div_x1, div_x2)) %*% c(div_x1, div_x2)) * (sum(diag(M_fun(motif, x1, x2))))/2)
    )
    ) %>%
  ungroup() %>% dplyr::select(motif, x1, x2, e_d, frequency)

ggplot(df_f, aes(x=as.factor(frequency), y=e_d)) + geom_boxplot()+facet_wrap(~motif)+ggtitle("M-matrix at startpoints")
ggplot(df_f, aes(x=as.factor(frequency), y=e_d)) + geom_violin()+facet_wrap(~motif)
ggplot(filter(df_f, x1 > 5, x2 > 5), aes(x=as.factor(frequency), y=e_d)) + geom_violin()+facet_wrap(~motif)
ggplot(filter(df_f, x1 > 5, x2 > 5), aes(x=as.factor(frequency), y=e_d)) + geom_boxplot()+facet_wrap(~motif)



# End points instead of starting points

df_f2 <- final_results %>%
  rowwise() %>%
  mutate(
    e_d = as.numeric(
      t(c(div_x1, div_x2)) %*% M_fun(motif, mean_X1, mean_X2) %*% c(div_x1, div_x2) /
        ((t(c(div_x1, div_x2)) %*% c(div_x1, div_x2)) * (sum(diag(M_fun(motif, mean_X1, mean_X2))))/2)
    )
  ) %>%
  ungroup() %>% dplyr::select(motif, x1, x2, e_d, frequency)

ggplot(df_f2, aes(x=as.factor(frequency), y=e_d)) + geom_boxplot()+facet_wrap(~motif, labeller = labeller(
  motif = function(x)
    paste("motif", x)))+ggtitle("")+
  labs(x="Landscape complexity", y = "Mutational variance in the direction of divergence")

ggsave("M-alignment.png", width = 8, height = 5, bg = "white")

ggplot(df_f2, aes(x=as.factor(frequency), y=e_d)) + geom_violin()+facet_wrap(~motif)

#################
# Angles

M <- M_fun(1, 5,5)
M
eig <- eigen(M1)

# Find index of largest eigenvalue (by magnitude)
idx <- which.max(abs(eig$values))

# Leading eigenvector
leading_vec <- eig$vectors[, idx]

leading_vec

eigen(M1)$vectors[, which.max(abs(eigen(M1)$values))]

df_angle <- final_results %>%
  rowwise() %>%
  mutate(
    angle = {
      vec1 = eigen(M_fun(motif, mean_X1, mean_X2) )$vectors[, which.max(abs(eigen(M_fun(motif, mean_X1, mean_X2) )$values))]
      
      acos(sum(vec1*c(div_x1, div_x2)) / (sqrt(sum(vec1^2)) * sqrt(sum(c(div_x1, div_x2)^2))))
    }
 
  ) %>%
  ungroup() %>% dplyr::select(motif, x1, x2, angle, frequency)

df_angle <- filter(final_results) %>%
  rowwise() %>%
  mutate(
    angle = {
      eig = eigen(M_fun(motif, mean_X1, mean_X2))
      vec1 = eig$vectors[, which.max(abs(eig$values))]
      
      v2 = c(div_x1, div_x2)
      
      cos_theta = sum(vec1 * v2) / (sqrt(sum(vec1^2)) * sqrt(sum(v2^2)))
      
      acos(abs(cos_theta))  # <-- key change
    },
    angle_norm = (angle / (pi / 2) ) * 90
  ) %>%
  ungroup() %>%
  dplyr::select(motif, x1, x2, angle, angle_norm, frequency, landscape)

ggplot(df_angle, aes(x=as.factor(frequency), y=angle_norm)) + geom_boxplot()+facet_wrap(~motif, labeller = labeller(
  motif = function(x)
    paste("motif", x)))+ggtitle("")+
  labs(x="Landscape complexity", y = "Angle between divergence and leading eigenvector")

ggplot(df_angle, aes(x = as.factor(motif),
                     y = angle_norm,
                     fill = as.factor(frequency))) +
  geom_boxplot(aes(group = interaction(motif, frequency)),
               position = position_dodge(width = 0.8)) +
  scale_y_continuous(
    breaks = c(0, 15, 30, 45, 60, 75, 90),
    labels = c("0°", "15°","30°","45°", "60°","75°", "90°")
  )+geom_hline(yintercept = 45, colour="red",linetype = "dashed")+
  labs(x = "Motif",
       y = "Angle between divergence and leading eigenvector",
       fill = "Frequency") +
  ggtitle("")+theme_bw()

ggplot(df_angle, aes(x = as.factor(frequency),
                     y = angle_norm,
                     fill = as.factor(motif))) +
  geom_boxplot(aes(group = interaction(frequency, motif)),
               position = position_dodge(width = 0.8)) +
  scale_y_continuous(
    breaks = c(0, 15, 30, 45, 60, 75, 90),
    labels = c("0°", "15°","30°","45°", "60°","75°", "90°")
  )+geom_hline(yintercept = 45, colour="red",linetype = "dashed")+
  labs(x = "Landscape complexity",
       y = expression("Angle between divergence and" ~ M[max]),
       fill = "Motif") +
  ggtitle("")+theme_bw()



ggsave("M-alignment_angle.png", width = 8, height = 5,  bg = "white", dpi = 600)





#######################################
# Example covariance matrix (2D)
M <- matrix(c(10, 9,
              9, 1), nrow=2, byrow=TRUE)
M <- M_fun(1, 10,5)

# Compute eigenvalues
eigs <- eigen(M)$values
cat("Eigenvalues:", eigs, "\n")

# Mean evolvability
k <- nrow(M)
e_mean <- sum(diag(M)) / k
cat("Mean evolvability:", e_mean, "\n")

# Function to compute normalized evolvability for unit vector at angle theta
normalized_evolvability <- function(theta, M, e_mean) {
  d <- c(cos(theta), sin(theta))
  e_d <- t(d) %*% M %*% d
  return(as.numeric(e_d / e_mean))
}

# Sequence of angles
angles <- seq(0, 2*pi, length.out = 360)
ratios <- sapply(angles, normalized_evolvability, M=M, e_mean=e_mean)

# Plot
plot(angles*180/pi, ratios, type='l', lwd=2,
     xlab="Angle (degrees)", ylab="Normalized evolvability (e_d / e_mean)",
     main="Normalized Evolvability vs Direction")
abline(h=1, col="red", lty=2)  # mean line






