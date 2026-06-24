

# Parameters for blended method
t1 <- 24
t2 <- 60
a <- 5
b <- 5

# Time points
time_horizon <- 84 
dt = 0.1
t_seq <- seq(0, time_horizon, dt)

plot_sens_beta <- function(t1, t2, a, b){
  # Compute blended/internal/external hazard
  h_Pem_selected <- extract_survival_hazard(m_Pem_selected, t_seq)$h
  h_Scha_selected <- extract_survival_hazard(m_Scha_selected, t_seq, rebased_time)$h
  h_Pem_blended <- compute_blended_hazard(h_Pem_selected, h_Scha_selected, t1, t2, a, b, t_seq) 
  
  # Plot blended hazard (vs fitted internal & external hazard)
  hazplot_Pem_blended <- plot_blended_hazard(h_Pem_blended, h_Pem_selected, h_Scha_selected, rebased_time, t1, t2, 'Pembrolizumab')
  hazplot_Pem_blended <- hazplot_Pem_blended + ggtitle(paste0("Beta(", a, ",", b, ")"))
  return(hazplot_Pem_blended)
}

betap1 <- plot_sens_beta(14,60,0.2,0.2)
betap2 <- plot_sens_beta(14,60,5,5)
betap3 <- plot_sens_beta(14,60,7,3)
betap4 <- plot_sens_beta(14,60,3,7)

betap1 + betap2 + betap3 + betap4 + 
  plot_layout(nrow = 2, byrow = TRUE, guides = "collect") & # "guides = 'collect'" puts all legends into one
  theme(legend.position = "bottom") # Position the collected legend




plot_sens_interval <- function(t1, t2, a, b){
  # Compute blended/internal/external hazard
  h_Pem_selected <- extract_survival_hazard(m_Pem_selected, t_seq)$h
  h_Scha_selected <- extract_survival_hazard(m_Scha_selected, t_seq, rebased_time)$h
  h_Pem_blended <- compute_blended_hazard(h_Pem_selected, h_Scha_selected, t1, t2, a, b, t_seq) 
  
  # Plot blended hazard (vs fitted internal & external hazard)
  hazplot_Pem_blended <- plot_blended_hazard(h_Pem_blended, h_Pem_selected, h_Scha_selected, rebased_time, t1, t2, 'Pembrolizumab')
  hazplot_Pem_blended <- hazplot_Pem_blended + ggtitle(paste0("Blending interval (", t1, ",", t2, ")"))
  return(hazplot_Pem_blended)
}

intervalp1 <- plot_sens_interval(14,36,3,7)
intervalp2 <- plot_sens_interval(14,60,3,7)
intervalp3 <- plot_sens_interval(24,36,3,7)
intervalp4 <- plot_sens_interval(24,60,3,7)

intervalp1 + intervalp2 + intervalp3 + intervalp4 + 
  plot_layout(nrow = 2, byrow = TRUE, guides = "collect") & # "guides = 'collect'" puts all legends into one
  theme(legend.position = "bottom") # Position the collected legend
